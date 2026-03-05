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

    auto isProc = [](const word& w)
    {
        return w.size() >= 9 && w.substr(0,9) == "processor";
    };

    auto isTime = [](const word& w)
    {
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
        {
            p = p.path();
            continue;
        }
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
    Wm_(0.0), Sm_(0.0),
    Wcv_(0.0), Scv_(0.0),
    bufferMean_(),
    bufferCV_(),
    outPath_(""),
    kpiPath_(""),
    lastInstMeanS_(0.0),
    lastInstStdS_(0.0),
    lastInstCVs_(0.0)
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

// compute area-weighted mean and std over patch (global via gSum)
void GamaPatchAvgTimeWindow::patchAreaMeanStd
(
    const volScalarField& f,
    scalar& mean,
    scalar& stdS
) const
{
    const fvPatch& p = this->mesh_.boundary()[patchi_];
    const scalarField& vals = f.boundaryField()[patchi_];
    const scalarField& Af   = p.magSf();

    const scalar A  = gSum(Af);
    const scalar m1 = (A > SMALL) ? gSum(vals*Af)/A : 0.0;
    const scalar m2 = (A > SMALL) ? gSum(vals*vals*Af)/A : 0.0;

    const scalar var = max(m2 - m1*m1, scalar(0));
    mean = m1;
    stdS = std::sqrt(var);
}

void GamaPatchAvgTimeWindow::pushSample
(
    std::deque<Sample>& buf,
    scalar& W,
    scalar& S,
    scalar dt,
    scalar value
)
{
    if (dt <= SMALL) return;
    buf.push_back({dt, value});
    W += dt;
    S += value * dt;
}

void GamaPatchAvgTimeWindow::trimWindow
(
    std::deque<Sample>& buf,
    scalar& W,
    scalar& S
)
{
    while (W - timeWindow_ > SMALL && !buf.empty())
    {
        scalar excess = W - timeWindow_;
        Sample& front = buf.front();

        if (front.dt <= excess + SMALL)
        {
            W -= front.dt;
            S -= front.value * front.dt;
            buf.pop_front();
        }
        else
        {
            const scalar cut = excess;
            front.dt -= cut;
            W -= cut;
            S -= front.value * cut;
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

    scalar meanS = 0.0;
    scalar stdS  = 0.0;
    patchAreaMeanStd(vf, meanS, stdS);

    lastInstMeanS_ = meanS;
    lastInstStdS_  = stdS;
    lastInstCVs_   = (mag(meanS) > SMALL) ? (stdS/meanS) : 0.0;

    scalar dt = mesh.time().deltaTValue();
    reduce(dt, minOp<scalar>());

    // moving window means
    pushSample(bufferMean_, Wm_,  Sm_,  dt, meanS);
    trimWindow(bufferMean_, Wm_,  Sm_);

    pushSample(bufferCV_,   Wcv_, Scv_, dt, lastInstCVs_);
    trimWindow(bufferCV_,   Wcv_, Scv_);

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
                    << "# Columns: Time\tmeanS\tstdS\tCVs\twinMean(meanS)\twinMean(CVs)" << nl;
            out_->flush();
        }
    }

    const scalar tNow = T.value();

    const scalar instMeanS = lastInstMeanS_;
    const scalar instStdS  = lastInstStdS_;
    const scalar instCVs   = lastInstCVs_;

    const bool haveWinM  = (Wm_  > SMALL);
    const bool haveWinCV = (Wcv_ > SMALL);

    const scalar winMeanM  = haveWinM  ? (Sm_/Wm_)   : instMeanS;
    const scalar winMeanCV = haveWinCV ? (Scv_/Wcv_) : instCVs;

    out_->setf(std::ios::scientific, std::ios::floatfield);
    (*out_) << tNow << '\t'
            << instMeanS << '\t' << instStdS << '\t' << instCVs << '\t'
            << winMeanM  << '\t' << winMeanCV << nl;
    out_->flush();

    // KPI = window-mean of CVs if available, else instantaneous CVs
    std::ofstream kpi(kpiPath_.c_str(), std::ios::out | std::ios::trunc);
    if (kpi.good())
    {
        kpi.setf(std::ios::scientific);
        kpi.precision(10);
        kpi << winMeanCV << '\n';
        kpi.flush();
    }
    else
    {
        Warning<< "Failed to open KPI file for writing: " << kpiPath_ << nl;
    }

    return true;
}

} // namespace Foam
