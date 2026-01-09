#include "membraneSolventFluxFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "OSspecific.H"
#include "OFstream.H"
#include "PstreamReduceOps.H"

#include <fstream>
#include <iomanip>

namespace Foam
{

defineTypeNameAndDebug(membraneSolventFluxFvPatchVectorField,0);
addToRunTimeSelectionTable
(
    fvPatchVectorField,
    membraneSolventFluxFvPatchVectorField,
    dictionary
);

// ------------------------------------------------------------------------- //
membraneSolventFluxFvPatchVectorField::
membraneSolventFluxFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector,volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p,iF)
{
    this->operator==(vector::zero);
    Jv_.setSize(p.size(),0.0);
}

membraneSolventFluxFvPatchVectorField::
membraneSolventFluxFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector,volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p,iF)
{
    sigma_   = dict.lookupOrDefault<scalar>("sigma",1.0);
    CAName_  = dict.lookupOrDefault<word>("CAName","CA");
    CApName_ = dict.lookupOrDefault<word>("CApName","CAp");

    word m(dict.lookupOrDefault<word>("model","SpieglerKedem"));
    if      (m=="SpieglerKedem")     model_=Model::SpieglerKedem;
    else if (m=="constantFlux")      model_=Model::constantFlux;
    else if (m=="targetAverageFlux") model_=Model::targetAverageFlux;
    else                             model_=Model::Unknown;

    if (model_==Model::SpieglerKedem) Ah_=dict.lookupOrDefault<scalar>("Ah",0.0);
    if (model_==Model::constantFlux)  JvConst_=dict.lookupOrDefault<scalar>("Jv",0.0);
    if (model_==Model::targetAverageFlux)
    {
        JvAverage_=dict.lookupOrDefault<scalar>("JvAverage",0.0);
        AhMin_    =dict.lookupOrDefault<scalar>("AhMin",0.0);
        AhMax_    =dict.lookupOrDefault<scalar>("AhMax",VGREAT);
    }

    // virialCoeffs: polynomial in Pa with C in kg/m^3 (no T-scaling)
    if (dict.found("virialCoeffs"))
    {
        virialCoeffs_ = List<scalar>(dict.lookup("virialCoeffs"));
    }
    if ((model_==Model::SpieglerKedem || model_==Model::targetAverageFlux)
        && virialCoeffs_.empty())
    {
        FatalIOErrorInFunction(dict)
            << "Missing 'virialCoeffs' for model "
            << (model_==Model::SpieglerKedem ? "SpieglerKedem" : "targetAverageFlux")
            << ". Provide 'virialCoeffs ( a1 a2 ... );' with C in kg/m^3 and π in Pa."
            << exit(FatalIOError);
    }

    // Δp always uses pPermConst_ (default 0)
    pPermConst_ = dict.lookupOrDefault<scalar>("pPermConst",0.0);

    // Optional moving window average
    windowAverage_ = dict.lookupOrDefault<Switch>("windowAverage",false);
    Twindow_       = dict.lookupOrDefault<scalar>("Twindow",0.0);
    if (!windowAverage_ || Twindow_ <= SMALL) { windowAverage_=false; Twindow_=0.0; }

    this->operator==(vector::zero);
    if (dict.found("value"))
        fvPatchVectorField::operator=(vectorField("value",dict,p.size()));
    Jv_.setSize(p.size(),0.0);
}

membraneSolventFluxFvPatchVectorField::
membraneSolventFluxFvPatchVectorField
(
    const membraneSolventFluxFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector,volMesh>& iF,
    const fvPatchFieldMapper& m
)
:
    fixedValueFvPatchVectorField(ptf,p,iF,m),
    pPtr_(nullptr), CAPtr_(nullptr), CApPtr_(nullptr),
    Ah_(ptf.Ah_), sigma_(ptf.sigma_), JvConst_(ptf.JvConst_),
    virialCoeffs_(ptf.virialCoeffs_),
    pPermConst_(ptf.pPermConst_),
    JvAverage_(ptf.JvAverage_), AhMin_(ptf.AhMin_), AhMax_(ptf.AhMax_),
    lastAhEff_(ptf.lastAhEff_), windowAverage_(ptf.windowAverage_),
    Twindow_(ptf.Twindow_), lastJvAvgInstant_(ptf.lastJvAvgInstant_),
    lastJvAvgWindow_(ptf.lastJvAvgWindow_), tBuf_(ptf.tBuf_), jvBuf_(ptf.jvBuf_),
    Jv_(ptf.Jv_), CAName_(ptf.CAName_), CApName_(ptf.CApName_), model_(ptf.model_)
{}

membraneSolventFluxFvPatchVectorField::
membraneSolventFluxFvPatchVectorField(const membraneSolventFluxFvPatchVectorField& ptf)
:
    fixedValueFvPatchVectorField(ptf),
    pPtr_(nullptr), CAPtr_(nullptr), CApPtr_(nullptr),
    Ah_(ptf.Ah_), sigma_(ptf.sigma_), JvConst_(ptf.JvConst_),
    virialCoeffs_(ptf.virialCoeffs_),
    pPermConst_(ptf.pPermConst_),
    JvAverage_(ptf.JvAverage_), AhMin_(ptf.AhMin_), AhMax_(ptf.AhMax_),
    lastAhEff_(ptf.lastAhEff_), windowAverage_(ptf.windowAverage_),
    Twindow_(ptf.Twindow_), lastJvAvgInstant_(ptf.lastJvAvgInstant_),
    lastJvAvgWindow_(ptf.lastJvAvgWindow_), tBuf_(ptf.tBuf_), jvBuf_(ptf.jvBuf_),
    Jv_(ptf.Jv_), CAName_(ptf.CAName_), CApName_(ptf.CApName_), model_(ptf.model_)
{}

membraneSolventFluxFvPatchVectorField::
membraneSolventFluxFvPatchVectorField
(
    const membraneSolventFluxFvPatchVectorField& ptf,
    const DimensionedField<vector,volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf,iF),
    pPtr_(nullptr), CAPtr_(nullptr), CApPtr_(nullptr),
    Ah_(ptf.Ah_), sigma_(ptf.sigma_), JvConst_(ptf.JvConst_),
    virialCoeffs_(ptf.virialCoeffs_),
    pPermConst_(ptf.pPermConst_),
    JvAverage_(ptf.JvAverage_), AhMin_(ptf.AhMin_), AhMax_(ptf.AhMax_),
    lastAhEff_(ptf.lastAhEff_), windowAverage_(ptf.windowAverage_),
    Twindow_(ptf.Twindow_), lastJvAvgInstant_(ptf.lastJvAvgInstant_),
    lastJvAvgWindow_(ptf.lastJvAvgWindow_), tBuf_(ptf.tBuf_), jvBuf_(ptf.jvBuf_),
    Jv_(ptf.Jv_), CAName_(ptf.CAName_), CApName_(ptf.CApName_), model_(ptf.model_)
{}

// ------------------------------------------------------------------------- //
scalar membraneSolventFluxFvPatchVectorField::
osmoticPressure_(scalar C) const
{
    // π(C) = Σ_k a_k C^k  with C in kg/m^3, π in Pa
    scalar sum = 0.0;
    scalar powC = C; // start at k=1
    for (label k=0; k<virialCoeffs_.size(); ++k)
    {
        sum += virialCoeffs_[k] * powC;
        powC *= C;
    }
    return sum;
}

void membraneSolventFluxFvPatchVectorField::
updateWindowAveraging_(const vectorField& Up) const
{
    if(!windowAverage_)
    {
        lastJvAvgWindow_ = lastJvAvgInstant_;
        return;
    }

    tmp<vectorField> tnf = patch().nf();
    const vectorField& n = tnf();
    const scalarField& Sf = patch().magSf();

    scalar num=0.0, den=0.0;
    forAll(Up,i){ const scalar Jv = (Up[i] & n[i]); num += Jv*Sf[i]; den += Sf[i]; }
    reduce(num,sumOp<scalar>()); reduce(den,sumOp<scalar>());
    lastJvAvgInstant_ = (den>SMALL ? num/den : 0.0);

    const scalar tNow = this->db().time().value();
    tBuf_.append(tNow); jvBuf_.append(lastJvAvgInstant_);

    if (Twindow_ <= SMALL)
    {
        const label maxKeep=1024;
        if (tBuf_.size() > maxKeep)
        { const label drop = tBuf_.size() - maxKeep;
          for (label k=0; k<drop; ++k) { tBuf_.remove(0); jvBuf_.remove(0); } }
        lastJvAvgWindow_ = lastJvAvgInstant_;
        return;
    }

    const scalar tMin = tNow - Twindow_;
    while (tBuf_.size()>0 && tBuf_[0] < tMin) { tBuf_.remove(0); jvBuf_.remove(0); }

    if (tBuf_.size() < 2) { lastJvAvgWindow_ = lastJvAvgInstant_; return; }

    const scalar tStart=tBuf_[0], tEnd=tBuf_[tBuf_.size()-1];
    scalar width = min(Twindow_, max(0.0, tEnd - tStart));
    if (width <= SMALL) { lastJvAvgWindow_ = lastJvAvgInstant_; return; }

    scalar integral=0.0;
    for (label k=1; k<tBuf_.size(); ++k)
    {
        const scalar dt = tBuf_[k]-tBuf_[k-1];
        const scalar avg= 0.5*(jvBuf_[k]+jvBuf_[k-1]);
        if (dt > 0) integral += dt*avg;
    }
    lastJvAvgWindow_ = integral/width;
}

// -------- NEW helper: append area-averaged permeate flux line (master-only) -
static void appendPermeateFluxLine_
(
    const word& patchName,
    const Time& runTime,
    const scalarField& Jv,
    const scalarField& Sf
)
{
    // area-average with parallel reduction
    scalar num = 0.0, den = 0.0;
    forAll(Jv, i) { num += Jv[i]*Sf[i]; den += Sf[i]; }
    reduce(num, sumOp<scalar>()); reduce(den, sumOp<scalar>());
    const scalar JvAvg = (den > SMALL ? num/den : 0.0);

    if (!Pstream::master()) return;

    const fileName dir = runTime.globalPath()/"postProcessing"/"membraneSolventFlux"/patchName;

    mkDir(dir);

    const fileName f = dir/"permeateFlux.dat";
    const bool existed = isFile(f);

    std::ofstream ofs(f.c_str(), std::ios_base::app);
    if (!ofs.good()) return;

    if (!existed)
    {
        ofs << "# Permeate flux at write times for patch " << patchName << "\n"
            << "# time[s]    Jv_avg[m/s]\n";
    }

    ofs.setf(std::ios::scientific);
    ofs.precision(10);
    ofs << std::setw(16) << runTime.value()
        << "  " << std::setw(16) << JvAvg << "\n";
}
// ------------------------------------------------------------------------- //

void membraneSolventFluxFvPatchVectorField::updateCoeffs()
{
    if (updated()) return;

    const label patchI = patch().index();

    // Look up required fields
    if (!pPtr_)   pPtr_   = &db().lookupObject<volScalarField>("p");
    if (!CAPtr_)  CAPtr_  = &db().lookupObject<volScalarField>(CAName_);
    if (!CApPtr_) CApPtr_ = &db().lookupObject<volScalarField>(CApName_);

    // Fail-fast: boundary sizes must match this patch
    auto checkSize = [&](const char* name, const volScalarField& f)
    {
        const label sz = f.boundaryField()[patchI].size();
        if (sz != patch().size())
        {
            FatalErrorInFunction
                << "Field '" << name << "' boundary size (" << sz
                << ") != patch size (" << patch().size() << ") on patch '"
                << patch().name() << "'. Ensure matching boundary entries."
                << exit(FatalError);
        }
    };
    checkSize("p",   *pPtr_);
    checkSize(CAName_.c_str(),  *CAPtr_);
    checkSize(CApName_.c_str(), *CApPtr_);

    tmp<vectorField> tnf = patch().nf();
    const vectorField& n = tnf();
    vectorField Up(patch().size(), vector::zero);

    switch (model_)
    {
        case Model::SpieglerKedem:
        {
            const scalarField& pf = pPtr_->boundaryField()[patchI];   // feed p
            const scalarField& Cf = CAPtr_->boundaryField()[patchI];  // feed CA
            const scalarField& Cp = CApPtr_->boundaryField()[patchI]; // perm CAp

            scalarField piF(Cf.size(), 0.0), piP(Cp.size(), 0.0);
            forAll(piF,i) { piF[i] = osmoticPressure_(Cf[i]); }
            forAll(piP,i) { piP[i] = osmoticPressure_(Cp[i]); }

            scalarField Jv(pf.size(), 0.0);
            forAll(Jv,i)
            {
                const scalar dp  = pf[i] - pPermConst_;   // Δp uses pPermConst_
                const scalar dpi = piF[i] - piP[i];       // Δπ uses CA/CAp
                Jv[i] = max(0.0, Ah_ * (dp - sigma_ * dpi)); // no retroflux
            }
            forAll(Up,i) Up[i] = n[i] * Jv[i];
            Jv_ = Jv;
            break;
        }

        case Model::constantFlux:
        {
            scalarField Jv(patch().size(), max(0.0, JvConst_)); // ensure ≥0
            forAll(Up,i) Up[i] = n[i] * Jv[i];
            Jv_ = Jv;
            break;
        }

        case Model::targetAverageFlux:
        {
            const scalarField& pf = pPtr_->boundaryField()[patchI];
            const scalarField& Cf = CAPtr_->boundaryField()[patchI];
            const scalarField& Cp = CApPtr_->boundaryField()[patchI];
            const scalarField& Sf = patch().magSf();

            scalarField drive(pf.size(), 0.0);
            forAll(drive,i)
            {
                const scalar piF = osmoticPressure_(Cf[i]);
                const scalar piP = osmoticPressure_(Cp[i]);
                const scalar dp  = pf[i] - pPermConst_;
                const scalar dpi = piF - piP;
                drive[i] = dp - sigma_ * dpi;
            }

            scalar num=0.0, den=0.0;
            forAll(drive,i){ num += Sf[i]*drive[i]; den += Sf[i]; }
            reduce(num,sumOp<scalar>()); reduce(den,sumOp<scalar>());
            const scalar meanDrive = (den>SMALL ? num/den : 0.0);

            scalar AhEff = 0.0;
            if (mag(meanDrive) > SMALL)
                AhEff = max(AhMin_, min(JvAverage_/meanDrive, AhMax_));
            else
                AhEff = 0.0;

            lastAhEff_ = AhEff;

            forAll(Up,i)
            {
                const scalar JvLocal = max(0.0, AhEff * drive[i]); // no retroflux
                Up[i] = n[i] * JvLocal;
                Jv_[i] = JvLocal;
            }
            break;
        }

        default:
            FatalErrorInFunction << "Unknown solvent flux model" << exit(FatalError);
    }

    this->operator==(Up);
    updateWindowAveraging_(Up);
    fixedValueFvPatchVectorField::updateCoeffs();
}

// ------------------------------------------------------------------------- //
void membraneSolventFluxFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);

    word m("Unknown");
    if      (model_==Model::SpieglerKedem)     m="SpieglerKedem";
    else if (model_==Model::constantFlux)      m="constantFlux";
    else if (model_==Model::targetAverageFlux) m="targetAverageFlux";

    os.writeKeyword("model")       << m        << token::END_STATEMENT << nl;
    os.writeKeyword("sigma")       << sigma_   << token::END_STATEMENT << nl;
    os.writeKeyword("CAName")      << CAName_  << token::END_STATEMENT << nl;
    os.writeKeyword("CApName")     << CApName_ << token::END_STATEMENT << nl;
    os.writeKeyword("virialCoeffs")<< virialCoeffs_ << token::END_STATEMENT << nl;

    // Δp reference (always written, default 0)
    os.writeKeyword("pPermConst")  << pPermConst_ << token::END_STATEMENT << nl;

    if (model_==Model::SpieglerKedem)
        os.writeKeyword("Ah") << Ah_ << token::END_STATEMENT << nl;
    else if (model_==Model::constantFlux)
        os.writeKeyword("Jv") << JvConst_ << token::END_STATEMENT << nl;
    else if (model_==Model::targetAverageFlux)
    {
        os.writeKeyword("JvAverage") << JvAverage_ << token::END_STATEMENT << nl;
        os.writeKeyword("AhMin")     << AhMin_     << token::END_STATEMENT << nl;
        os.writeKeyword("AhMax")     << AhMax_     << token::END_STATEMENT << nl;
        os.writeKeyword("AhStar")    << lastAhEff_ << token::END_STATEMENT << nl;
    }

    os.writeKeyword("windowAverage") << Switch(windowAverage_) << token::END_STATEMENT << nl;
    if (windowAverage_)
    {
        os.writeKeyword("Twindow")       << Twindow_           << token::END_STATEMENT << nl;
        os.writeKeyword("JvAvg_instant") << lastJvAvgInstant_ << token::END_STATEMENT << nl;
        os.writeKeyword("JvAvg_window")  << lastJvAvgWindow_  << token::END_STATEMENT << nl;
    }

    writeEntry("value", os);

    // NEW: write permeate flux log at write times (SpieglerKedem only)
    if (model_ == Model::SpieglerKedem)
    {
        const Time& runTime = this->db().time();
        const scalarField& Sf = patch().magSf();
        // Jv_ holds outward-positive per-face flux [m/s]
        appendPermeateFluxLine_(patch().name(), runTime, Jv_, Sf);
    }

    // Master-only summary for targetAverageFlux (existing)
if (model_==Model::targetAverageFlux && Pstream::master())
{
    const Time& runTime = this->db().time();

    // ANTES:
    // fileName dir = runTime.path()/"postProcessing"/"membraneSolventFlux"/patch().name();

    // DEPOIS (raiz do case, fora de processor*/):
    fileName dir = runTime.globalPath()/"postProcessing"/"membraneSolventFlux"/patch().name();

    mkDir(dir);
    OFstream osCalib(dir/"calibratedProperties.dat");
    if (osCalib.good())
    {
        osCalib << "# Calibrated membrane properties for patch " << patch().name() << nl
                << "timeFinal     " << runTime.timeOutputValue() << nl
                << "AhStar        " << lastAhEff_                << nl
                << "JvAverage     " << JvAverage_                << nl
                << "sigma         " << sigma_                    << nl
                << "virialCoeffs  " << virialCoeffs_             << nl
                << "pPermConst    " << pPermConst_               << nl
                << "model         targetAverageFlux"             << nl;
    }
}

}

} // namespace Foam
// ************************************************************************* //
