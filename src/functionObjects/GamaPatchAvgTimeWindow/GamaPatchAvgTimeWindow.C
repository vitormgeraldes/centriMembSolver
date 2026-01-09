#include "GamaPatchAvgTimeWindow.H"
#include "addToRunTimeSelectionTable.H"
#include "IOstreams.H"
#include "OSspecific.H"
#include "Pstream.H"
#include "PstreamReduceOps.H"
#include "Time.H"
#include "volFields.H"

#include <deque>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <sys/stat.h>
#include <cctype>

namespace Foam
{

defineTypeNameAndDebug(GamaPatchAvgTimeWindow, 0);
addToRunTimeSelectionTable(functionObject, GamaPatchAvgTimeWindow, dictionary);

// ---------------- utility: file size ----------------
static bool fileExistsSize(const fileName& f, off_t& sz)
{
    struct stat sb;
    if (::stat(f.c_str(), &sb) == 0)
    {
        sz = sb.st_size;
        return true;
    }
    sz = 0;
    return false;
}

// ---------------- utility: find case root (strip processor*/time) ----------------
static inline fileName caseRoot(const Time& T)
{
    fileName p = T.rootPath()/T.caseName();

    auto isProc = [](const word& w){ return w.size() >= 9 && w.substr(0,9) == "processor"; };
    auto isTime = [](const word& w){
        if (w.empty()) return false;
        bool hasDigit = false;
        for (char c : std::string(w))
        {
            if (std::isdigit(static_cast<unsigned char>(c))) hasDigit = true;
            else if (c != '.' && c != '-' && c != 'e' && c != 'E' && c != '+') return false;
        }
        return hasDigit;
    };

    for (;;)
    {
        const word tail = p.name();
        if (isProc(tail) || isTime(tail) || tail == "constant" || tail == "system")
        { p = p.path(); continue; }
        break;
    }
    return p;
}

// ============================================================================

GamaPatchAvgTimeWindow::GamaPatchAvgTimeWindow
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::fvMeshFunctionObject(name, runTime, dict),
    patchName_(dict.lookup("patchName")),
    fieldName_(dict.lookup("fieldName")),
    outputName_(dict.lookupOrDefault<word>("outputName","Gama_membrane_avgWin")),
    fileName_(dict.lookupOrDefault<word>("fileName","areaAvg_timeWindow.dat")),
    timeWindow_(readScalar(dict.lookup("timeWindow"))),
    masterOnly_(dict.lookupOrDefault<Switch>("masterOnly", true)),
    patchi_(-1),
    W_(0.0), S1_(0.0), S2_(0.0),
    outPath_(""),
    kpiPath_(""),
    wroteHeader_(false),
    lastInstAvg_(0.0)
{
    if (timeWindow_ <= SMALL)
    {
        FatalErrorInFunction
            << "timeWindow must be > 0 (s). Got " << timeWindow_ << nl
            << exit(FatalError);
    }

    const fvMesh& mesh = this->mesh_;
    patchi_ = mesh.boundaryMesh().findPatchID(patchName_);
    if (patchi_ < 0)
    {
        FatalErrorInFunction
            << "Patch '" << patchName_ << "' not found. Available: "
            << mesh.boundaryMesh().names() << nl << exit(FatalError);
    }

    const fileName root = caseRoot(mesh.time());

    if (Pstream::master())
    {
        if (!isDir(root/"postProcessing"))             mkDir(root/"postProcessing");
        if (!isDir(root/"postProcessing"/outputName_)) mkDir(root/"postProcessing"/outputName_);
        if (!isDir(root/"postProcessing"/"KPI"))       mkDir(root/"postProcessing"/"KPI");
    }

    outPath_ = root/"postProcessing"/outputName_/fileName_;
    kpiPath_ = root/"postProcessing"/"KPI"/("Gama_" + patchName_ + ".txt");

    if (Pstream::master())
    {
        Info<< "GamaPatchAvgTimeWindow: outPath=" << outPath_
            << " | kpiPath=" << kpiPath_ << nl;
    }
}

// compute area-weighted average over patch (global via gSum)
scalar GamaPatchAvgTimeWindow::patchAreaAverage(const volScalarField& f) const
{
    const fvPatch& p = this->mesh_.boundary()[patchi_];
    const scalarField& vals = f.boundaryField()[patchi_];
    const scalarField& Af = p.magSf();

    const scalar num = gSum(vals * Af);
    const scalar den = gSum(Af);
    return den > SMALL ? num/den : 0.0;
}

void GamaPatchAvgTimeWindow::pushSample(scalar dt, scalar value)
{
    if (dt <= SMALL) return;
    buffer_.push_back({dt, value});
    W_  += dt;
    S1_ += value * dt;
    S2_ += value * value * dt;
}

void GamaPatchAvgTimeWindow::trimWindow()
{
    while (W_ - timeWindow_ > SMALL && !buffer_.empty())
    {
        scalar excess = W_ - timeWindow_;
        Sample& front = buffer_.front();
        if (front.dt <= excess + SMALL)
        {
            W_  -= front.dt;
            S1_ -= front.value * front.dt;
            S2_ -= front.value * front.value * front.dt;
            buffer_.pop_front();
        }
        else
        {
            const scalar cut = excess;
            front.dt -= cut;
            W_  -= cut;
            S1_ -= front.value * cut;
            S2_ -= front.value * front.value * cut;
        }
    }
}

// ---------------- execute: all ranks participate ----------------
bool GamaPatchAvgTimeWindow::execute()
{
    const fvMesh& mesh = this->mesh_;

    bool have = mesh.foundObject<volScalarField>(fieldName_);
    reduce(have, andOp<bool>());
    if (!have) return true;

    const volScalarField& vf = mesh.lookupObject<volScalarField>(fieldName_);
    const scalar areaAvg = patchAreaAverage(vf);
    lastInstAvg_ = areaAvg;    // store for master-only write

    scalar dt = mesh.time().deltaTValue();
    reduce(dt, minOp<scalar>());

    pushSample(dt, areaAvg);
    trimWindow();
    return true;
}

// ---------------- write: master-only, no MPI collectives ----------------
bool GamaPatchAvgTimeWindow::write()
{
    if (masterOnly_ && !Pstream::master()) return true;

    const Time& T = this->mesh_.time();
    const fileName root = caseRoot(T);

    if (Pstream::master())
    {
        if (!isDir(root/"postProcessing"))             mkDir(root/"postProcessing");
        if (!isDir(root/"postProcessing"/outputName_)) mkDir(root/"postProcessing"/outputName_);
        if (!isDir(root/"postProcessing"/"KPI"))       mkDir(root/"postProcessing"/"KPI");
    }

    off_t sz = 0;
    bool existed = fileExistsSize(outPath_, sz);
    if (!out_.valid())
    {
        IOstreamOption opt(IOstreamOption::ASCII, IOstreamOption::UNCOMPRESSED);
        out_.reset(new OFstream(outPath_, opt,
            existed ? IOstreamOption::APPEND : IOstreamOption::NO_APPEND));

        if (!out_->good())
        {
            Warning<< "Failed to open output file: " << outPath_ << nl;
        }
        out_->precision(10);
        if (!existed || sz == 0)
        {
            (*out_) << "# Patch: " << patchName_ << nl
                    << "# Field: " << fieldName_ << nl
                    << "# timeWindow(s): " << timeWindow_ << nl
                    << "# Columns: Time\tareaAverage(" << fieldName_ << "@"
                    << patchName_ << ")\twinMean\twinVar\twinStd" << nl;
            out_->flush();
        }
    }

    const scalar inst = lastInstAvg_;
    const scalar tNow = T.value();

    scalar mean=0, var=0, stdv=0;
    const bool haveWin = (W_ > SMALL);
    if (haveWin)
    {
        mean = S1_/W_;
        const scalar m2 = S2_/W_ - mean*mean;
        var  = m2 > 0 ? m2 : 0.0;
        stdv = std::sqrt(var);
    }

    out_->setf(std::ios::scientific, std::ios::floatfield);
    (*out_) << tNow << '\t' << inst << '\t';
    if (haveWin) (*out_) << mean << '\t' << var << '\t' << stdv << nl;
    else         (*out_) << "n/a\tn/a\tn/a" << nl;
    out_->flush();

    // write KPI (overwrite)
    std::ofstream kpi(kpiPath_.c_str(), std::ios::out | std::ios::trunc);
    if (kpi.good())
    {
        kpi.setf(std::ios::scientific);
        kpi.precision(10);
        kpi << (haveWin ? mean : inst) << '\n';
        kpi.flush();
    }
    else
    {
        Warning<< "Failed to open KPI file for writing: " << kpiPath_ << nl;
    }

    return true;
}

} // namespace Foam
