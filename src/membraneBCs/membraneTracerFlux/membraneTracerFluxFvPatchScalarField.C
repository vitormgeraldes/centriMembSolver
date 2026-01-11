/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Website:  www.openfoam.com                      |
|   \\  /    A nd           | Version:  v2506                                 |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
/*  membraneTracerFluxFvPatchScalarField
 *
 *  Tracer-flux membrane boundary condition (scalar field).
 *  Couples diffusive transport with convection at a selective interface.
 *
 *  Models implemented (Spiegler–Kedem / Kedem–Katchalsky family):
 *   1) intrinsicRejection
 *        JsT = Jv * CTp,  with CTp = (1 - RTint) * CTw
 *        Film model (Robin BC on CA): k(CTw - CTi) = (RTint * Jv) * CTw
 *        → CTw = [kT / (kT - RTint*Jv)] * CTi
 *
 *   2) tracerPermeability  (SK closure)
 *        SK equation: JsT = (1 - sigmaT) Jv CTw + BT (CTw - CTp)
 *        → CTp = [(BT + (1 - sigmaT)Jv)/(BT + Jv)] * CTw = αT*CTw
 *        and JsT = Jv * CTp
 *        If 'sigmaT' omitted, σ is inherited from the solvent BC on U (defaults to 1).
 *
 *   3) observedRejection (inverse)
 *        Infers RTint so that RTobs_model = 1 - <Jv*CTp>/<Jv>/CTb matches user 'RTobs'
 *        Uses a safeguarded secant+bisection. MPI-safe reductions.
 *
 *  Diagnostics (if CTb > 0):
 *   GamaT = (CTw - CTb)/CTb
 *   RTobserved = 1 - <Jv * CTp> / <Jv> / CTb
 *
 *  Output:
 *   - Writes per-face JsT, Jv, CTp, GamaT to the patch dictionary
 *   - Logs time-series in postProcessing/membraneTracerFlux/<patch>/...
 *   - Summary tables in postProcessing/membraneTracerFlux/_summary
 *
 *  Units:
 *   CT, CTb [kg/m3];  Jv [m/s];  JsT [kg/m2/s];  DTeff [m2/s];
 *   BT [m/s];  sigmaT, RTint, RTobs, GamaT [-]
 *
 *  Notes:
 *   - Parallel/periodic safe (uses returnReduce, master-only logging guards)
 *   - Stable for small denominators (kT - RTint*Jv), clamps when needed
 *   - For stationary meshes (no mesh motion assumed)
 *   - This BC may optionally “mirror” its per-face diagnostics into existing
 *     volScalarField(s) named "CTp", "JsT", "Jv" if present (see end of updateCoeffs()).
 */
/*---------------------------------------------------------------------------*/

#include "membraneTracerFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "membraneSolventFluxFvPatchVectorField.H"
#include "volFields.H"
#include "PstreamReduceOps.H"
#include "Pstream.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "Time.H"

#include <fstream>   // std::ofstream
#include <iomanip>   // std::setw, std::setprecision

namespace Foam
{

// Runtime registration
defineTypeNameAndDebug(membraneTracerFluxFvPatchScalarField, 0);
addToRunTimeSelectionTable
(
    fvPatchScalarField,
    membraneTracerFluxFvPatchScalarField,
    dictionary
);

// -------------------- Constructors --------------------

membraneTracerFluxFvPatchScalarField::membraneTracerFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    JsT_(p.size(), Zero),
    Jv_(p.size(), Zero),
    CTp_(p.size(), Zero),
    GamaT_(p.size(), Zero),
    RTobserved_(0.0),
    lastLoggedTimeIndex_(-1)
{}

membraneTracerFluxFvPatchScalarField::membraneTracerFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    JsT_(p.size(), Zero),
    Jv_(p.size(), Zero),
    CTp_(p.size(), Zero),
    GamaT_(p.size(), Zero),
    RTobserved_(0.0),
    lastLoggedTimeIndex_(-1)
{
    // Field name for velocity (to fetch solvent BC and normals·U)
    UName_ = dict.lookupOrDefault<word>("UName", "U");

    // Select model
    const word m(dict.lookupOrDefault<word>("model", "intrinsicRejection"));
    if      (m == "intrinsicRejection") model_ = Model::intrinsicRejection;
    else if (m == "solutePermeability") model_ = Model::solutePermeability;
    else if (m == "observedRejection")  model_ = Model::observedRejection;
    else                                model_ = Model::Unknown;

    // Reference/bulk concentration for diagnostics and inverse mode
    CTb_ = dict.lookupOrDefault<scalar>("CTb", 0.0);

    // Model parameters
    if (model_ == Model::intrinsicRejection)
    {
        RTint_ = dict.lookupOrDefault<scalar>("RTint", 0.0);
    }
    else if (model_ == Model::solutePermeability)
    {
        BT_ = dict.lookupOrDefault<scalar>("BT", 0.0);
        hasSigmaT_ = dict.found("sigmaT");
        if (hasSigmaT_) sigmaT_ = readScalar(dict.lookup("sigmaT"));
    }
    else if (model_ == Model::observedRejection)
    {
        RTobs_    = dict.lookupOrDefault<scalar>("RTobs", 0.0);
        RTintMin_ = dict.lookupOrDefault<scalar>("RTintMin", 0.0);
        RTintMax_ = dict.lookupOrDefault<scalar>("RTintMax", 0.9999);
        tol_     = dict.lookupOrDefault<scalar>("tol", 1e-6);
        maxIter_ = dict.lookupOrDefault<label>("maxIter", 50);

        // Minimal seeding if Rint not specified: start from Robs
        if (!dict.found("RTint")) RTint_ = RTobs_;
    }

    // Optional initial "value" (standard for patch fields)
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
    }
}

// Copy-like constructors
membraneTracerFluxFvPatchScalarField::membraneTracerFluxFvPatchScalarField
(
    const membraneTracerFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& m
)
:
    mixedFvPatchScalarField(ptf, p, iF, m),
    UName_(ptf.UName_),
    RTint_(ptf.RTint_),
    BT_(ptf.BT_),
    sigmaT_(ptf.sigmaT_),
    hasSigmaT_(ptf.hasSigmaT_),
    RTobs_(ptf.RTobs_),
    RTobsAchieved_(ptf.RTobsAchieved_),
    CTb_(ptf.CTb_),
    RTintMin_(ptf.RTintMin_),
    RTintMax_(ptf.RTintMax_),
    tol_(ptf.tol_),
    maxIter_(ptf.maxIter_),
    model_(ptf.model_),
    JsT_(ptf.JsT_),
    Jv_(ptf.Jv_),
    CTp_(ptf.CTp_),
    GamaT_(ptf.GamaT_),
    RTobserved_(ptf.RTobserved_),
    lastLoggedTimeIndex_(ptf.lastLoggedTimeIndex_)
{}

membraneTracerFluxFvPatchScalarField::membraneTracerFluxFvPatchScalarField
(
    const membraneTracerFluxFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    UName_(ptf.UName_),
    RTint_(ptf.RTint_),
    BT_(ptf.BT_),
    sigmaT_(ptf.sigmaT_),
    hasSigmaT_(ptf.hasSigmaT_),
    RTobs_(ptf.RTobs_),
    RTobsAchieved_(ptf.RTobsAchieved_),
    CTb_(ptf.CTb_),
    RTintMin_(ptf.RTintMin_),
    RTintMax_(ptf.RTintMax_),
    tol_(ptf.tol_),
    maxIter_(ptf.maxIter_),
    model_(ptf.model_),
    JsT_(ptf.JsT_),
    Jv_(ptf.Jv_),
    CTp_(ptf.CTp_),
    GamaT_(ptf.GamaT_),
    RTobserved_(ptf.RTobserved_),
    lastLoggedTimeIndex_(ptf.lastLoggedTimeIndex_)
{}

membraneTracerFluxFvPatchScalarField::membraneTracerFluxFvPatchScalarField
(
    const membraneTracerFluxFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    UName_(ptf.UName_),
    RTint_(ptf.RTint_),
    BT_(ptf.BT_),
    sigmaT_(ptf.sigmaT_),
    hasSigmaT_(ptf.hasSigmaT_),
    RTobs_(ptf.RTobs_),
    RTobsAchieved_(ptf.RTobsAchieved_),
    CTb_(ptf.CTb_),
    RTintMin_(ptf.RTintMin_),
    RTintMax_(ptf.RTintMax_),
    tol_(ptf.tol_),
    maxIter_(ptf.maxIter_),
    model_(ptf.model_),
    JsT_(ptf.JsT_),
    Jv_(ptf.Jv_),
    CTp_(ptf.CTp_),
    GamaT_(ptf.GamaT_),
    RTobserved_(ptf.RTobserved_),
    lastLoggedTimeIndex_(ptf.lastLoggedTimeIndex_)
{}

// -------------------- Helpers --------------------

scalar membraneTracerFluxFvPatchScalarField::sigmaTFromUOrDefault1_
(
    const label patchI
) const
{
    // Try to read σ from the membraneSolventFlux BC on U (same patch).
    // Fall back to 1.0 if not available.
    const volVectorField& U = db().lookupObject<volVectorField>(UName_);
    const fvPatchVectorField& Ubc = U.boundaryField()[patchI];

    if (const auto* solv =
            dynamic_cast<const membraneSolventFluxFvPatchVectorField*>(&Ubc))
    {
        return solv->sigma();
    }
    return 1.0;
}

scalar membraneTracerFluxFvPatchScalarField::evaluateForRTint_
(
    const scalar RTint,
    const scalarField& CTi,
    const scalarField& kTDel,
    const scalarField& Jv,
    const scalarField& Ap,
    scalar& CTpFluxAvg
) const
{
    // Compute flux-weighted CTp average for a given Rint (LOCAL only here).
    const label n = CTi.size();
    scalar sumJvA = 0.0, sumJvCTpA = 0.0;

    for (label i=0; i<n; ++i)
    {
        const scalar k = max(kTDel[i], SMALL);
        const scalar J = Jv[i];
        const scalar denom = k - RTint*J;
        const scalar CTw = (denom > SMALL ? (k/denom)*CTi[i] : CTi[i]);
        const scalar CTp = max(0.0, (1.0 - RTint) * CTw);

        sumJvA    += J * Ap[i];
        sumJvCTpA  += J * CTp * Ap[i];
    }

    CTpFluxAvg = (sumJvA > VSMALL ? sumJvCTpA/sumJvA : 0.0);
    return CTpFluxAvg;
}

scalar membraneTracerFluxFvPatchScalarField::inferRTintFromObserved_
(
    const scalarField& CTi,
    const scalarField& kTDel,
    const scalarField& Jv,
    const scalarField& Ap,
    scalar& RTobsAchieved
) const
{
    // RTobs(model) at R, computed from GLOBAL (reduced) sums each call
    auto RTobsOf = [&](const scalar RT)->scalar
    {
        scalar sumJvA_local   = 0.0;
        scalar sumJvCTpA_local = 0.0;
        const label n = CTi.size();

        for (label i=0; i<n; ++i)
        {
            const scalar k = max(kTDel[i], SMALL);
            const scalar J = Jv[i];
            const scalar denom = k - RT*J;
            const scalar CTw = (denom > SMALL ? (k/denom)*CTi[i] : CTi[i]);
            const scalar CTp = max(0.0, (1.0 - RT) * CTw);

            sumJvA_local   += J * Ap[i];
            sumJvCTpA_local += J * CTp * Ap[i];
        }

        const scalar sumJvA   = returnReduce(sumJvA_local,   sumOp<scalar>());
        const scalar sumJvCTpA = returnReduce(sumJvCTpA_local, sumOp<scalar>());

        const scalar CTpFluxAvg_global =
            (sumJvA > VSMALL ? sumJvCTpA / sumJvA : 0.0);

        return 1.0 - (CTb_ > SMALL ? CTpFluxAvg_global / CTb_ : 0.0);
    };

    const scalar target = RTobs_;
    const scalar aG = clamp(RTintMin_, 0.0, 0.999999);
    const scalar bG = clamp(RTintMax_, aG + SMALL, 0.9999999);

    // Seed from last step
    scalar RT = clamp(RTint_, aG, bG);
    scalar f = RTobsOf(RT);
    if (mag(f - target) < tol_) { RTobsAchieved = f; return RT; }

    // Safeguarded secant parameters
    scalar d0  = 0.01*(bG - aG);   // probe step (1%)
    const scalar dMax = 0.10*(bG - aG);

    for (int it=0; it<maxIter_; ++it)
    {
        // Symmetric probe to estimate slope
        scalar d = min(d0, min(RT - aG, bG - RT) * 0.5);
        if (d <= SMALL) d = 0.5*(bG - aG)*1e-3;

        const scalar RTm = max(aG, RT - d);
        const scalar RTp = min(bG, RT + d);
        const scalar fm = RTobsOf(RTm);
        const scalar fp = RTobsOf(RTp);
        const scalar slope = (fp - fm) / max(RTp - RTm, SMALL);

        if (mag(slope) < VSMALL)
        {
            // Too flat: gentle nudge towards target
            const scalar dir  = (f < target ? 1.0 : -1.0);
            const scalar step = clamp(0.25*dMax, 1e-5, dMax);
            const scalar RTtry = clamp(RT + dir*step, aG, bG);
            const scalar ftry = RTobsOf(RTtry);

            if (mag(ftry - target) < mag(f - target))
            {
                RT = RTtry; f = ftry;
                if (mag(f - target) < tol_) break;
            }
            else
            {
                // Robust fallback: global bisection
                scalar a = aG, b = bG, fa = RTobsOf(a), fb = RTobsOf(b);

                if ((fa - target) > 0 && (fb - target) > 0)
                { RTobsAchieved = fb; return b; }
                if ((fa - target) < 0 && (fb - target) < 0)
                { RTobsAchieved = fa; return a; }

                for (int it2=0; it2<maxIter_; ++it2)
                {
                    const scalar mid = 0.5*(a + b);
                    const scalar fm2 = RTobsOf(mid);
                    if (mag(fm2 - target) < tol_) { RTobsAchieved = fm2; return mid; }

                    if ((fa - target)*(fm2 - target) <= 0)
                    { b = mid; fb = fm2; }
                    else
                    { a = mid; fa = fm2; }
                }
                const scalar mid = 0.5*(a + b);
                RTobsAchieved = RTobsOf(mid);
                return mid;
            }
            continue;
        }

        // Secant/Newton-like update (clamped)
        scalar RTnext = RT + (target - f)/slope;
        const scalar jump = clamp(RTnext - RT, -dMax, dMax);
        RTnext = clamp(RT + jump, aG, bG);

        const scalar fnext = RTobsOf(RTnext);
        if (mag(fnext - target) < tol_) { RTobsAchieved = fnext; return RTnext; }

        if (mag(fnext - target) < mag(f - target))
        { RT = RTnext; f = fnext; }
        else
        { d0 *= 0.5; }
    }

    RTobsAchieved = f;
    return RT;
}

// -------------------- updateCoeffs --------------------

void membraneTracerFluxFvPatchScalarField::updateCoeffs()
{
    if (updated()) return;

    const label patchI = patch().index();
// --- FAIL-FAST: exigir campos de exportação (CTp, Js, Jv) ---
{
    const label pid = patchI;

    auto requireField = [&](const char* name, const char* dims) -> const volScalarField&
    {
        // 1) Tem de existir no objectRegistry
        if (!db().foundObject<volScalarField>(name))
        {
            FatalErrorInFunction
                << "Required volScalarField '" << name << "' not found for patch '"
                << patch().name() << "'.\n"
                << "Create it (e.g., in createFields.H or 0/) with dimensions " << dims << ".\n"
                << "Example:\n"
                << "  volScalarField " << name << "(IOobject(\"" << name << "\", runTime.timeName(), mesh,\n"
                << "      IOobject::NO_READ, IOobject::AUTO_WRITE), mesh,\n"
                << "      dimensionedScalar(\"" << name << "\", " << dims << ", 0.0));\n"
                << exit(FatalError);
        }

        // 2) E o tamanho do boundary no patch tem de coincidir
        const volScalarField& f = db().lookupObject<volScalarField>(name);
        if (f.boundaryField()[pid].size() != patch().size())
        {
            FatalErrorInFunction
                << "Field '" << name << "' boundary size ("
                << f.boundaryField()[pid].size() << ") != patch size ("
                << patch().size() << ") on patch '" << patch().name() << "'.\n"
                << "Ensure the field has a boundary entry on this patch with matching size."
                << exit(FatalError);
        }

        return f;
    };

    //   CTp : [1 -3 0 0 0 0 0]  (kg/m3)
    //   JsT : [1 -2 -1 0 0 0 0] (kg/m2/s)
    //   Jv  : [0  1 -1 0 0 0 0] (m/s)
    (void)requireField("CTp", "dimensionSet(1,-3,0,0,0,0,0)");
    (void)requireField("JsT",  "dimensionSet(1,-2,-1,0,0,0,0)");
    (void)requireField("Jv",  "dimensionSet(0, 1,-1,0,0,0,0)");
}
// --- fim FAIL-FAST ---

    // Geometry & fields
    const scalarField& delta = patch().deltaCoeffs();    // [1/m]
    const scalarField& Ap    = patch().magSf();          // [m2]
    const vectorField  n     = patch().nf();

    const volVectorField& U = db().lookupObject<volVectorField>(UName_);
    const vectorField& Up   = U.boundaryField()[patchI];

    const volScalarField& DTeff = db().lookupObject<volScalarField>("DTeff");
    const scalarField& DTp      = DTeff.boundaryField()[patchI]; // [m2/s]

    const scalarField Ci = this->patchInternalField();   // [kg/m3] (owner side)

    // Outward-positive volumetric flux Jv [m/s]
    Jv_.setSize(patch().size(), 0.0);
    forAll(Jv_, i) Jv_[i] = max(Up[i] & n[i], scalar(0));

    // Film conductance k = D * delta [m/s]
    scalarField kTDel(patch().size(), 0.0);
    forAll(kTDel, i) kTDel[i] = max(DTp[i]*delta[i], SMALL);

    // Prepare model state
    scalar RTintUsed = RTint_;
    RTobsAchieved_ = 0.0;

    if (model_ == Model::observedRejection)
    {
        // Require some permeation and a positive CAb before inverse
        scalar localSumJvA = 0.0;
        forAll(Jv_, i) localSumJvA += Jv_[i] * Ap[i];
        const scalar sumJvA = returnReduce(localSumJvA, sumOp<scalar>());

        if (sumJvA <= VSMALL || CTb_ <= SMALL)
        {
            WarningInFunction
                << "observedRejection: insufficient permeation or CTb<=0 on patch '"
                << patch().name() << "'. Keeping RTint=" << RTint_ << nl;
            RTintUsed = RTint_;
        }
        else
        {
            RTintUsed = inferRTintFromObserved_(Ci, kTDel, Jv_, Ap, RTobsAchieved_);
        }
    }

    // Mixed (pure Robin) initialization
    this->refValue()      = scalarField(patch().size(), 0.0);
    this->refGrad()       = scalarField(patch().size(), 0.0);
    this->valueFraction() = scalarField(patch().size(), 0.0);

    // Per-face diagnostics
    JsT_.setSize(patch().size(), 0.0);
    CTp_.setSize(patch().size(), 0.0);
    GamaT_.setSize(patch().size(), 0.0);

    // sigmaT used (only for solutePermeability)
    const scalar sigmaTEffDefault =
        (model_ == Model::solutePermeability)
        ? (hasSigmaT_ ? sigmaT_ : sigmaTFromUOrDefault1_(patchI))
        : 1.0;

    // Face loop
    forAll(Ci, i)
    {
        const scalar k = max(kTDel[i], SMALL);   // [m/s]
        const scalar J = Jv_[i];                // [m/s]

        // hMix is the "Robin" term we push to the LHS (negative of physical coupling)
        scalar hMix = 0.0;
        if (model_ == Model::intrinsicRejection || model_ == Model::observedRejection)
        {
            const scalar hPhys = RTintUsed * J;
            hMix = -hPhys;
        }
        else if (model_ == Model::solutePermeability)
        {
            const scalar denomBTJ = max(BT_ + J, SMALL);
            const scalar hPhys   = (J > SMALL ? sigmaTEffDefault * J*J / denomBTJ : 0.0);
            hMix = -hPhys;
        }
        else
        {
            FatalErrorInFunction
                << "Unknown model on patch " << patch().name() << nl
                << "Valid: intrinsicRejection | solutePermeability | observedRejection"
                << exit(FatalError);
        }

        // Build the mixed BC coefficients
        const scalar denom = k + hMix; // = k - hPhys
        scalar vf = 0.0;
        scalar CTw = Ci[i];

        if (denom > SMALL)
        {
            vf = hMix/denom;        // valueFraction
            CTw = (k/denom)*Ci[i];   // wall concentration
        }
        else
        {
            // Degenerate case: keep numerically safe
            const scalar denomSafe = SMALL;
            vf = hMix/denomSafe;
            CTw = (k/denomSafe)*Ci[i];
        }
        this->valueFraction()[i] = vf;

        // Permeate concentration Cp and solute flux Js
        if (model_ == Model::intrinsicRejection || model_ == Model::observedRejection)
        {
            const scalar CTp = max(0.0, (1.0 - RTintUsed) * CTw);
            CTp_[i] = CTp;
            JsT_[i]  = J * CTp_[i];
        }
        else // solutePermeability
        {
            const scalar alpha =
                (BT_ + (1.0 - sigmaTEffDefault)*J) / max(BT_ + J, SMALL);
            const scalar CTp = max(0.0, alpha * CTw);
            CTp_[i] = CTp;
            JsT_[i]  = J * CTp_[i];
        }

        // Polarization Γ = (CTw - CAb)/CAb (if CAb>0)
        GamaT_[i] = (CTb_ > SMALL) ? (CTw - CTb_)/CTb_ : 0.0;
    }

    // Observed rejection from flux-weighted mean (global, MPI-safe)
    if (CTb_ > SMALL)
    {
        scalar localSumJvA   = 0.0;
        scalar localSumJvCTpA = 0.0;
        forAll(Jv_, i)
        {
            localSumJvA   += Jv_[i] * Ap[i];
            localSumJvCTpA += Jv_[i]*CTp_[i]* Ap[i];
        }
        const scalar sumJvA   = returnReduce(localSumJvA,   sumOp<scalar>());
        const scalar sumJvCTpA = returnReduce(localSumJvCTpA, sumOp<scalar>());

        if (sumJvA > VSMALL)
        {
            const scalar CTpFluxAvg_global = sumJvCTpA/sumJvA;
            RTobserved_ = 1.0 - CTpFluxAvg_global / CTb_;
            RTobserved_ = max(0.0, min(1.0, RTobserved_));
        }
        else
        {
            RTobserved_ = 0.0; // negligible permeation
        }
    }
    else
    {
        RTobserved_ = 0.0; // undefined without reference → keep numeric
    }

    // Persist the RTint used by the inverse model (so it is written out)
    if (model_ == Model::observedRejection)
    {
        const_cast<membraneTracerFluxFvPatchScalarField*>(this)->RTint_ = RTintUsed;
    }

    // (Optional) Export diagnostics directly into existing volScalarFields.
    // The BC does not create fields; it only fills them if present.
    try
    {
        const label pid = patch().index();

        if (db().foundObject<volScalarField>("CTp"))
        {
            auto& CTpField =
                const_cast<volScalarField&>(db().lookupObject<volScalarField>("CTp"));
            if (CTpField.boundaryField()[pid].size() == CTp_.size())
            {
                CTpField.boundaryFieldRef()[pid] = CTp_;
            }
        }
        if (db().foundObject<volScalarField>("JsT"))
        {
            auto& JsTField =
                const_cast<volScalarField&>(db().lookupObject<volScalarField>("JsT"));
            if (JsTField.boundaryField()[pid].size() == JsT_.size())
            {
                JsTField.boundaryFieldRef()[pid] = JsT_;
            }
        }
        if (db().foundObject<volScalarField>("Jv"))
        {
            auto& JvField =
                const_cast<volScalarField&>(db().lookupObject<volScalarField>("Jv"));
            if (JvField.boundaryField()[pid].size() == Jv_.size())
            {
                JvField.boundaryFieldRef()[pid] = Jv_;
            }
        }
    }
    catch(const std::exception& e)
    {
        WarningInFunction << "Export to volScalarField(s) failed: " << e.what() << nl;
    }

    mixedFvPatchScalarField::updateCoeffs();
}

// -------------------- write --------------------

void membraneTracerFluxFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    // Model name for output
    word modelName("Unknown");
    switch (model_)
    {
        case Model::intrinsicRejection:  modelName = "intrinsicRejection";  break;
        case Model::solutePermeability:  modelName = "solutePermeability";  break;
        case Model::observedRejection:   modelName = "observedRejection";   break;
        default: break;
    }

    os.writeKeyword("model")  << modelName << token::END_STATEMENT << nl;
    os.writeKeyword("UName")  << UName_    << token::END_STATEMENT << nl;

    if (model_ == Model::intrinsicRejection || model_ == Model::observedRejection)
    {
        os.writeKeyword("RTint") << RTint_ << token::END_STATEMENT << nl;
    }

    if (model_ == Model::solutePermeability)
    {
        const scalar sigmaTUsed =
            (hasSigmaT_ ? sigmaT_ : sigmaTFromUOrDefault1_(patch().index()));
        os.writeKeyword("BT")          << BT_         << token::END_STATEMENT << nl;
        os.writeKeyword("sigmaTUsed")  << sigmaTUsed  << token::END_STATEMENT << nl;
        if (hasSigmaT_)
        {
            os.writeKeyword("sigmaTExplicit") << sigmaT_ << token::END_STATEMENT << nl;
        }
    }

    if (model_ == Model::observedRejection)
    {
        os.writeKeyword("RTobs")          << RTobs_          << token::END_STATEMENT << nl;
        os.writeKeyword("RTobsAchieved")  << RTobsAchieved_  << token::END_STATEMENT << nl;
    }

    // Reference and per-face diagnostics
    os.writeKeyword("CTb") << CTb_ << token::END_STATEMENT << nl;
    os.writeKeyword("JsT")  << JsT_  << token::END_STATEMENT << nl;
    os.writeKeyword("Jv")  << Jv_  << token::END_STATEMENT << nl;
    os.writeKeyword("CTp") << CTp_ << token::END_STATEMENT << nl;

    if (CTb_ > SMALL)
    {
        os.writeKeyword("GamaT")      << GamaT_      << token::END_STATEMENT << nl;
        os.writeKeyword("RTobserved") << RTobserved_ << token::END_STATEMENT << nl;
    }

    writeEntry("value", os);

    // Per-write logging (master only)
    if (Pstream::master())
    {
        const Time& runTime = this->db().time();
        if (runTime.outputTime() && runTime.timeIndex() != lastLoggedTimeIndex_)
        {
            lastLoggedTimeIndex_ = runTime.timeIndex();

            // Helper for case-root postProcessing dirs (serial & parallel)
            auto casePostDir = [&](const word& sub) -> fileName
            {
                return runTime.globalPath() / "postProcessing" / sub;
            };

            // Ensure base directories exist
            const fileName base    = casePostDir("membraneTracerFlux") / patch().name();
            const fileName sumBase = casePostDir("membraneTracerFlux") / "_summary";
            mkDir(base);
            mkDir(sumBase);

            // Per-patch time-series
            if (model_ == Model::solutePermeability)
            {
                const fileName fp = base/"BT_sigmaT_RTobserved_vs_time.dat";
                const bool existed = isFile(fp);
                std::ofstream osf(fp.c_str(), std::ios::out | std::ios::app);
                if (osf.good())
                {
                    if (!existed) osf << "# time[s]\tBT\tsigmaTUsed\tRTobserved\n";
                    const scalar sigmaTUsed =
                        (hasSigmaT_ ? sigmaT_ : sigmaTFromUOrDefault1_(patch().index()));
                    osf << runTime.value() << '\t' << BT_ << '\t'
                        << sigmaTUsed << '\t' << RTobserved_ << '\n';
                }
            }
            else
            {
                const fileName fp = base/"RTint_vs_time.dat";
                const bool existed = isFile(fp);
                std::ofstream osf(fp.c_str(), std::ios::out | std::ios::app);
                if (osf.good())
                {
                    if (!existed) osf << "# time[s]\tRTint\tRTobserved\n";
                    osf << runTime.value() << '\t' << RTint_ << '\t' << RTobserved_ << '\n';
                }
            }

            // Summary tables (_summary), overwritten each writeTime
            static label lastSummaryTimeIndexGlobal = -1;
            const bool newWriteTime =
                (runTime.timeIndex() != lastSummaryTimeIndexGlobal);

            if (model_ == Model::solutePermeability)
            {
                const fileName latest = sumBase/"BT_sigmaT_latest.dat";
                if (newWriteTime)
                {
                    std::ofstream osLatest(latest.c_str(), std::ios::out | std::ios::trunc);
                    if (osLatest.good())
                        osLatest << "# time[s]\tpatch\tBT\tsigmaTUsed\tRTobserved\n";
                    lastSummaryTimeIndexGlobal = runTime.timeIndex();
                }
                std::ofstream osLatest(latest.c_str(), std::ios::out | std::ios::app);
                if (osLatest.good())
                {
                    const scalar sigmaTUsed =
                        (hasSigmaT_ ? sigmaT_ : sigmaTFromUOrDefault1_(patch().index()));
                    osLatest << runTime.value() << '\t' << patch().name() << '\t'
                             << BT_ << '\t' << sigmaTUsed << '\t' << RTobserved_ << '\n';
                }

                const fileName snap = sumBase/("BT_sigmaT_" + runTime.timeName() + ".dat");
                const bool existedSnap = isFile(snap);
                std::ofstream osSnap(snap.c_str(), std::ios::out | std::ios::app);
                if (osSnap.good())
                {
                    if (!existedSnap)
                        osSnap << "# time[s]\tpatch\tBT\tsigmaTUsed\tRTobserved\n";
                    const scalar sigmaTUsed =
                        (hasSigmaT_ ? sigmaT_ : sigmaTFromUOrDefault1_(patch().index()));
                    osSnap << runTime.value() << '\t' << patch().name() << '\t'
                           << BT_ << '\t' << sigmaTUsed << '\t' << RTobserved_ << '\n';
                }
            }
            else
            {
                const fileName latest = sumBase/"RTint_latest.dat";
                if (newWriteTime)
                {
                    std::ofstream osLatest(latest.c_str(), std::ios::out | std::ios::trunc);
                    if (osLatest.good())
                        osLatest << "# time[s]\tpatch\tRTint\tRTobserved\n";
                    lastSummaryTimeIndexGlobal = runTime.timeIndex();
                }
                std::ofstream osLatest(latest.c_str(), std::ios::out | std::ios::app);
                if (osLatest.good())
                {
                    osLatest << runTime.value() << '\t' << patch().name() << '\t'
                             << RTint_ << '\t' << RTobserved_ << '\n';
                }

                const fileName snap = sumBase/("RTint_" + runTime.timeName() + ".dat");
                const bool existedSnap = isFile(snap);
                std::ofstream osSnap(snap.c_str(), std::ios::out | std::ios::app);
                if (osSnap.good())
                {
                    if (!existedSnap)
                        osSnap << "# time[s]\tpatch\tRTint\tRTobserved\n";
                    osSnap << runTime.value() << '\t' << patch().name() << '\t'
                           << RTint_ << '\t' << RTobserved_ << '\n';
                }
            }
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
