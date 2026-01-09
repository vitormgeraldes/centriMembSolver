/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Website:  www.openfoam.com                      |
|   \\  /    A nd           | Version:  v2506                                 |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
/*  membraneSoluteFluxFvPatchScalarField
 *
 *  Solute-flux membrane boundary condition (scalar field).
 *  Couples diffusive transport with convection at a selective interface.
 *
 *  Models implemented (Spiegler–Kedem / Kedem–Katchalsky family):
 *   1) intrinsicRejection
 *        Js = Jv * Cp,  with Cp = (1 - Rint) * Cw
 *        Film model (Robin BC on CA): k(Cw - Ci) = (Rint * Jv) * Cw
 *        → Cw = [k / (k - Rint*Jv)] * Ci
 *
 *   2) solutePermeability  (SK closure)
 *        SK equation: Js = (1 - σ) Jv Cw + B (Cw - Cp)
 *        → Cp = [(B + (1 - σ)Jv)/(B + Jv)] * Cw = α*Cw
 *        and Js = Jv * Cp
 *        If 'sigma' omitted, σ is inherited from the solvent BC on U (defaults to 1).
 *
 *   3) observedRejection (inverse)
 *        Infers Rint so that Robs_model = 1 - <Jv*Cp>/<Jv>/CAb matches user 'Robs'
 *        Uses a safeguarded secant+bisection. MPI-safe reductions.
 *
 *  Diagnostics (if CAb > 0):
 *   Γ = (Cw - CAb)/CAb
 *   Robserved = 1 - <Jv * Cp> / <Jv> / CAb
 *
 *  Output:
 *   - Writes per-face Js, Jv, CAp, Γ to the patch dictionary
 *   - Logs time-series in postProcessing/membraneSoluteFlux/<patch>/...
 *   - Summary tables in postProcessing/membraneSoluteFlux/_summary
 *
 *  Units:
 *   CA, CAb [kg/m3];  Jv [m/s];  Js [kg/m2/s];  Deff [m2/s];
 *   B [m/s];  σ, Rint, Robs, Γ [-]
 *
 *  Notes:
 *   - Parallel/periodic safe (uses returnReduce, master-only logging guards)
 *   - Stable for small denominators (k - Rint*Jv), clamps when needed
 *   - For stationary meshes (no mesh motion assumed)
 *   - This BC may optionally “mirror” its per-face diagnostics into existing
 *     volScalarField(s) named "CAp", "Js", "Jv" if present (see end of updateCoeffs()).
 */
/*---------------------------------------------------------------------------*/

#include "membraneSoluteFluxFvPatchScalarField.H"
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
defineTypeNameAndDebug(membraneSoluteFluxFvPatchScalarField, 0);
addToRunTimeSelectionTable
(
    fvPatchScalarField,
    membraneSoluteFluxFvPatchScalarField,
    dictionary
);

// -------------------- Constructors --------------------

membraneSoluteFluxFvPatchScalarField::membraneSoluteFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    Js_(p.size(), Zero),
    Jv_(p.size(), Zero),
    CAp_(p.size(), Zero),
    Gama_(p.size(), Zero),
    Robserved_(0.0),
    lastLoggedTimeIndex_(-1)
{}

membraneSoluteFluxFvPatchScalarField::membraneSoluteFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    Js_(p.size(), Zero),
    Jv_(p.size(), Zero),
    CAp_(p.size(), Zero),
    Gama_(p.size(), Zero),
    Robserved_(0.0),
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
    CAb_ = dict.lookupOrDefault<scalar>("CAb", 0.0);

    // Model parameters
    if (model_ == Model::intrinsicRejection)
    {
        Rint_ = dict.lookupOrDefault<scalar>("Rint", 0.0);
    }
    else if (model_ == Model::solutePermeability)
    {
        B_ = dict.lookupOrDefault<scalar>("B", 0.0);
        hasSigma_ = dict.found("sigma");
        if (hasSigma_) sigma_ = readScalar(dict.lookup("sigma"));
    }
    else if (model_ == Model::observedRejection)
    {
        Robs_    = dict.lookupOrDefault<scalar>("Robs", 0.0);
        RintMin_ = dict.lookupOrDefault<scalar>("RintMin", 0.0);
        RintMax_ = dict.lookupOrDefault<scalar>("RintMax", 0.9999);
        tol_     = dict.lookupOrDefault<scalar>("tol", 1e-6);
        maxIter_ = dict.lookupOrDefault<label>("maxIter", 50);

        // Minimal seeding if Rint not specified: start from Robs
        if (!dict.found("Rint")) Rint_ = Robs_;
    }

    // Optional initial "value" (standard for patch fields)
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
    }
}

// Copy-like constructors
membraneSoluteFluxFvPatchScalarField::membraneSoluteFluxFvPatchScalarField
(
    const membraneSoluteFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& m
)
:
    mixedFvPatchScalarField(ptf, p, iF, m),
    UName_(ptf.UName_),
    Rint_(ptf.Rint_),
    B_(ptf.B_),
    sigma_(ptf.sigma_),
    hasSigma_(ptf.hasSigma_),
    Robs_(ptf.Robs_),
    RobsAchieved_(ptf.RobsAchieved_),
    CAb_(ptf.CAb_),
    RintMin_(ptf.RintMin_),
    RintMax_(ptf.RintMax_),
    tol_(ptf.tol_),
    maxIter_(ptf.maxIter_),
    model_(ptf.model_),
    Js_(ptf.Js_),
    Jv_(ptf.Jv_),
    CAp_(ptf.CAp_),
    Gama_(ptf.Gama_),
    Robserved_(ptf.Robserved_),
    lastLoggedTimeIndex_(ptf.lastLoggedTimeIndex_)
{}

membraneSoluteFluxFvPatchScalarField::membraneSoluteFluxFvPatchScalarField
(
    const membraneSoluteFluxFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    UName_(ptf.UName_),
    Rint_(ptf.Rint_),
    B_(ptf.B_),
    sigma_(ptf.sigma_),
    hasSigma_(ptf.hasSigma_),
    Robs_(ptf.Robs_),
    RobsAchieved_(ptf.RobsAchieved_),
    CAb_(ptf.CAb_),
    RintMin_(ptf.RintMin_),
    RintMax_(ptf.RintMax_),
    tol_(ptf.tol_),
    maxIter_(ptf.maxIter_),
    model_(ptf.model_),
    Js_(ptf.Js_),
    Jv_(ptf.Jv_),
    CAp_(ptf.CAp_),
    Gama_(ptf.Gama_),
    Robserved_(ptf.Robserved_),
    lastLoggedTimeIndex_(ptf.lastLoggedTimeIndex_)
{}

membraneSoluteFluxFvPatchScalarField::membraneSoluteFluxFvPatchScalarField
(
    const membraneSoluteFluxFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    UName_(ptf.UName_),
    Rint_(ptf.Rint_),
    B_(ptf.B_),
    sigma_(ptf.sigma_),
    hasSigma_(ptf.hasSigma_),
    Robs_(ptf.Robs_),
    RobsAchieved_(ptf.RobsAchieved_),
    CAb_(ptf.CAb_),
    RintMin_(ptf.RintMin_),
    RintMax_(ptf.RintMax_),
    tol_(ptf.tol_),
    maxIter_(ptf.maxIter_),
    model_(ptf.model_),
    Js_(ptf.Js_),
    Jv_(ptf.Jv_),
    CAp_(ptf.CAp_),
    Gama_(ptf.Gama_),
    Robserved_(ptf.Robserved_),
    lastLoggedTimeIndex_(ptf.lastLoggedTimeIndex_)
{}

// -------------------- Helpers --------------------

scalar membraneSoluteFluxFvPatchScalarField::sigmaFromUOrDefault1_
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

scalar membraneSoluteFluxFvPatchScalarField::evaluateForRint_
(
    const scalar Rint,
    const scalarField& Ci,
    const scalarField& kDel,
    const scalarField& Jv,
    const scalarField& Ap,
    scalar& CpFluxAvg
) const
{
    // Compute flux-weighted Cp average for a given Rint (LOCAL only here).
    const label n = Ci.size();
    scalar sumJvA = 0.0, sumJvCpA = 0.0;

    for (label i=0; i<n; ++i)
    {
        const scalar k = max(kDel[i], SMALL);
        const scalar J = Jv[i];
        const scalar denom = k - Rint*J;
        const scalar Cw = (denom > SMALL ? (k/denom)*Ci[i] : Ci[i]);
        const scalar Cp = max(0.0, (1.0 - Rint) * Cw);

        sumJvA    += J * Ap[i];
        sumJvCpA  += J * Cp * Ap[i];
    }

    CpFluxAvg = (sumJvA > VSMALL ? sumJvCpA/sumJvA : 0.0);
    return CpFluxAvg;
}

scalar membraneSoluteFluxFvPatchScalarField::inferRintFromObserved_
(
    const scalarField& Ci,
    const scalarField& kDel,
    const scalarField& Jv,
    const scalarField& Ap,
    scalar& RobsAchieved
) const
{
    // Robs(model) at R, computed from GLOBAL (reduced) sums each call
    auto RobsOf = [&](const scalar R)->scalar
    {
        scalar sumJvA_local   = 0.0;
        scalar sumJvCpA_local = 0.0;
        const label n = Ci.size();

        for (label i=0; i<n; ++i)
        {
            const scalar k = max(kDel[i], SMALL);
            const scalar J = Jv[i];
            const scalar denom = k - R*J;
            const scalar Cw = (denom > SMALL ? (k/denom)*Ci[i] : Ci[i]);
            const scalar Cp = max(0.0, (1.0 - R) * Cw);

            sumJvA_local   += J * Ap[i];
            sumJvCpA_local += J * Cp * Ap[i];
        }

        const scalar sumJvA   = returnReduce(sumJvA_local,   sumOp<scalar>());
        const scalar sumJvCpA = returnReduce(sumJvCpA_local, sumOp<scalar>());

        const scalar CpFluxAvg_global =
            (sumJvA > VSMALL ? sumJvCpA / sumJvA : 0.0);

        return 1.0 - (CAb_ > SMALL ? CpFluxAvg_global / CAb_ : 0.0);
    };

    const scalar target = Robs_;
    const scalar aG = clamp(RintMin_, 0.0, 0.999999);
    const scalar bG = clamp(RintMax_, aG + SMALL, 0.9999999);

    // Seed from last step
    scalar R = clamp(Rint_, aG, bG);
    scalar f = RobsOf(R);
    if (mag(f - target) < tol_) { RobsAchieved = f; return R; }

    // Safeguarded secant parameters
    scalar d0  = 0.01*(bG - aG);   // probe step (1%)
    const scalar dMax = 0.10*(bG - aG);

    for (int it=0; it<maxIter_; ++it)
    {
        // Symmetric probe to estimate slope
        scalar d = min(d0, min(R - aG, bG - R) * 0.5);
        if (d <= SMALL) d = 0.5*(bG - aG)*1e-3;

        const scalar Rm = max(aG, R - d);
        const scalar Rp = min(bG, R + d);
        const scalar fm = RobsOf(Rm);
        const scalar fp = RobsOf(Rp);
        const scalar slope = (fp - fm) / max(Rp - Rm, SMALL);

        if (mag(slope) < VSMALL)
        {
            // Too flat: gentle nudge towards target
            const scalar dir  = (f < target ? 1.0 : -1.0);
            const scalar step = clamp(0.25*dMax, 1e-5, dMax);
            const scalar Rtry = clamp(R + dir*step, aG, bG);
            const scalar ftry = RobsOf(Rtry);

            if (mag(ftry - target) < mag(f - target))
            {
                R = Rtry; f = ftry;
                if (mag(f - target) < tol_) break;
            }
            else
            {
                // Robust fallback: global bisection
                scalar a = aG, b = bG, fa = RobsOf(a), fb = RobsOf(b);

                if ((fa - target) > 0 && (fb - target) > 0)
                { RobsAchieved = fb; return b; }
                if ((fa - target) < 0 && (fb - target) < 0)
                { RobsAchieved = fa; return a; }

                for (int it2=0; it2<maxIter_; ++it2)
                {
                    const scalar mid = 0.5*(a + b);
                    const scalar fm2 = RobsOf(mid);
                    if (mag(fm2 - target) < tol_) { RobsAchieved = fm2; return mid; }

                    if ((fa - target)*(fm2 - target) <= 0)
                    { b = mid; fb = fm2; }
                    else
                    { a = mid; fa = fm2; }
                }
                const scalar mid = 0.5*(a + b);
                RobsAchieved = RobsOf(mid);
                return mid;
            }
            continue;
        }

        // Secant/Newton-like update (clamped)
        scalar Rnext = R + (target - f)/slope;
        const scalar jump = clamp(Rnext - R, -dMax, dMax);
        Rnext = clamp(R + jump, aG, bG);

        const scalar fnext = RobsOf(Rnext);
        if (mag(fnext - target) < tol_) { RobsAchieved = fnext; return Rnext; }

        if (mag(fnext - target) < mag(f - target))
        { R = Rnext; f = fnext; }
        else
        { d0 *= 0.5; }
    }

    RobsAchieved = f;
    return R;
}

// -------------------- updateCoeffs --------------------

void membraneSoluteFluxFvPatchScalarField::updateCoeffs()
{
    if (updated()) return;

    const label patchI = patch().index();
// --- FAIL-FAST: exigir campos de exportação (CAp, Js, Jv) ---
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

    // Dimensões esperadas:
    //   CAp : [1 -3 0 0 0 0 0]  (kg/m3)
    //   Js  : [1 -2 -1 0 0 0 0] (kg/m2/s)
    //   Jv  : [0  1 -1 0 0 0 0] (m/s)
    (void)requireField("CAp", "dimensionSet(1,-3,0,0,0,0,0)");
    (void)requireField("Js",  "dimensionSet(1,-2,-1,0,0,0,0)");
    (void)requireField("Jv",  "dimensionSet(0, 1,-1,0,0,0,0)");
}
// --- fim FAIL-FAST ---

    // Geometry & fields
    const scalarField& delta = patch().deltaCoeffs();    // [1/m]
    const scalarField& Ap    = patch().magSf();          // [m2]
    const vectorField  n     = patch().nf();

    const volVectorField& U = db().lookupObject<volVectorField>(UName_);
    const vectorField& Up   = U.boundaryField()[patchI];

    const volScalarField& Deff = db().lookupObject<volScalarField>("Deff");
    const scalarField& Dp      = Deff.boundaryField()[patchI]; // [m2/s]

    const scalarField Ci = this->patchInternalField();   // [kg/m3] (owner side)

    // Outward-positive volumetric flux Jv [m/s]
    Jv_.setSize(patch().size(), 0.0);
    forAll(Jv_, i) Jv_[i] = max(Up[i] & n[i], scalar(0));

    // Film conductance k = D * delta [m/s]
    scalarField kDel(patch().size(), 0.0);
    forAll(kDel, i) kDel[i] = max(Dp[i]*delta[i], SMALL);

    // Prepare model state
    scalar RintUsed = Rint_;
    RobsAchieved_ = 0.0;

    if (model_ == Model::observedRejection)
    {
        // Require some permeation and a positive CAb before inverse
        scalar localSumJvA = 0.0;
        forAll(Jv_, i) localSumJvA += Jv_[i] * Ap[i];
        const scalar sumJvA = returnReduce(localSumJvA, sumOp<scalar>());

        if (sumJvA <= VSMALL || CAb_ <= SMALL)
        {
            WarningInFunction
                << "observedRejection: insufficient permeation or CAb<=0 on patch '"
                << patch().name() << "'. Keeping Rint=" << Rint_ << nl;
            RintUsed = Rint_;
        }
        else
        {
            RintUsed = inferRintFromObserved_(Ci, kDel, Jv_, Ap, RobsAchieved_);
        }
    }

    // Mixed (pure Robin) initialization
    this->refValue()      = scalarField(patch().size(), 0.0);
    this->refGrad()       = scalarField(patch().size(), 0.0);
    this->valueFraction() = scalarField(patch().size(), 0.0);

    // Per-face diagnostics
    Js_.setSize(patch().size(), 0.0);
    CAp_.setSize(patch().size(), 0.0);
    Gama_.setSize(patch().size(), 0.0);

    // σ used (only for solutePermeability)
    const scalar sigmaEffDefault =
        (model_ == Model::solutePermeability)
        ? (hasSigma_ ? sigma_ : sigmaFromUOrDefault1_(patchI))
        : 1.0;

    // Face loop
    forAll(Ci, i)
    {
        const scalar k = max(kDel[i], SMALL);   // [m/s]
        const scalar J = Jv_[i];                // [m/s]

        // hMix is the "Robin" term we push to the LHS (negative of physical coupling)
        scalar hMix = 0.0;
        if (model_ == Model::intrinsicRejection || model_ == Model::observedRejection)
        {
            const scalar hPhys = RintUsed * J;
            hMix = -hPhys;
        }
        else if (model_ == Model::solutePermeability)
        {
            const scalar denomBJ = max(B_ + J, SMALL);
            const scalar hPhys   = (J > SMALL ? sigmaEffDefault * J*J / denomBJ : 0.0);
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
        scalar Cw = Ci[i];

        if (denom > SMALL)
        {
            vf = hMix/denom;        // valueFraction
            Cw = (k/denom)*Ci[i];   // wall concentration
        }
        else
        {
            // Degenerate case: keep numerically safe
            const scalar denomSafe = SMALL;
            vf = hMix/denomSafe;
            Cw = (k/denomSafe)*Ci[i];
        }
        this->valueFraction()[i] = vf;

        // Permeate concentration Cp and solute flux Js
        if (model_ == Model::intrinsicRejection || model_ == Model::observedRejection)
        {
            const scalar Cp = max(0.0, (1.0 - RintUsed) * Cw);
            CAp_[i] = Cp;
            Js_[i]  = J * CAp_[i];
        }
        else // solutePermeability
        {
            const scalar alpha =
                (B_ + (1.0 - sigmaEffDefault)*J) / max(B_ + J, SMALL);
            const scalar Cp = max(0.0, alpha * Cw);
            CAp_[i] = Cp;
            Js_[i]  = J * CAp_[i];
        }

        // Polarization Γ = (Cw - CAb)/CAb (if CAb>0)
        Gama_[i] = (CAb_ > SMALL) ? (Cw - CAb_)/CAb_ : 0.0;
    }

    // Observed rejection from flux-weighted mean (global, MPI-safe)
    if (CAb_ > SMALL)
    {
        scalar localSumJvA   = 0.0;
        scalar localSumJvCpA = 0.0;
        forAll(Jv_, i)
        {
            localSumJvA   += Jv_[i] * Ap[i];
            localSumJvCpA += Jv_[i]*CAp_[i]* Ap[i];
        }
        const scalar sumJvA   = returnReduce(localSumJvA,   sumOp<scalar>());
        const scalar sumJvCpA = returnReduce(localSumJvCpA, sumOp<scalar>());

        if (sumJvA > VSMALL)
        {
            const scalar CpFluxAvg_global = sumJvCpA/sumJvA;
            Robserved_ = 1.0 - CpFluxAvg_global / CAb_;
            Robserved_ = max(0.0, min(1.0, Robserved_));
        }
        else
        {
            Robserved_ = 0.0; // negligible permeation
        }
    }
    else
    {
        Robserved_ = 0.0; // undefined without reference → keep numeric
    }

    // Persist the Rint used by the inverse model (so it is written out)
    if (model_ == Model::observedRejection)
    {
        const_cast<membraneSoluteFluxFvPatchScalarField*>(this)->Rint_ = RintUsed;
    }

    // (Optional) Export diagnostics directly into existing volScalarFields.
    // The BC does not create fields; it only fills them if present.
    try
    {
        const label pid = patch().index();

        if (db().foundObject<volScalarField>("CAp"))
        {
            auto& CApField =
                const_cast<volScalarField&>(db().lookupObject<volScalarField>("CAp"));
            if (CApField.boundaryField()[pid].size() == CAp_.size())
            {
                CApField.boundaryFieldRef()[pid] = CAp_;
            }
        }
        if (db().foundObject<volScalarField>("Js"))
        {
            auto& JsField =
                const_cast<volScalarField&>(db().lookupObject<volScalarField>("Js"));
            if (JsField.boundaryField()[pid].size() == Js_.size())
            {
                JsField.boundaryFieldRef()[pid] = Js_;
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

void membraneSoluteFluxFvPatchScalarField::write(Ostream& os) const
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
        os.writeKeyword("Rint") << Rint_ << token::END_STATEMENT << nl;
    }

    if (model_ == Model::solutePermeability)
    {
        const scalar sigmaUsed =
            (hasSigma_ ? sigma_ : sigmaFromUOrDefault1_(patch().index()));
        os.writeKeyword("B")          << B_         << token::END_STATEMENT << nl;
        os.writeKeyword("sigmaUsed")  << sigmaUsed  << token::END_STATEMENT << nl;
        if (hasSigma_)
        {
            os.writeKeyword("sigmaExplicit") << sigma_ << token::END_STATEMENT << nl;
        }
    }

    if (model_ == Model::observedRejection)
    {
        os.writeKeyword("Robs")          << Robs_          << token::END_STATEMENT << nl;
        os.writeKeyword("RobsAchieved")  << RobsAchieved_  << token::END_STATEMENT << nl;
    }

    // Reference and per-face diagnostics
    os.writeKeyword("CAb") << CAb_ << token::END_STATEMENT << nl;
    os.writeKeyword("Js")  << Js_  << token::END_STATEMENT << nl;
    os.writeKeyword("Jv")  << Jv_  << token::END_STATEMENT << nl;
    os.writeKeyword("CAp") << CAp_ << token::END_STATEMENT << nl;

    if (CAb_ > SMALL)
    {
        os.writeKeyword("Gama")      << Gama_      << token::END_STATEMENT << nl;
        os.writeKeyword("Robserved") << Robserved_ << token::END_STATEMENT << nl;
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
            const fileName base    = casePostDir("membraneSoluteFlux") / patch().name();
            const fileName sumBase = casePostDir("membraneSoluteFlux") / "_summary";
            mkDir(base);
            mkDir(sumBase);

            // Per-patch time-series
            if (model_ == Model::solutePermeability)
            {
                const fileName fp = base/"B_sigma_Robserved_vs_time.dat";
                const bool existed = isFile(fp);
                std::ofstream osf(fp.c_str(), std::ios::out | std::ios::app);
                if (osf.good())
                {
                    if (!existed) osf << "# time[s]\tB\tsigmaUsed\tRobserved\n";
                    const scalar sigmaUsed =
                        (hasSigma_ ? sigma_ : sigmaFromUOrDefault1_(patch().index()));
                    osf << runTime.value() << '\t' << B_ << '\t'
                        << sigmaUsed << '\t' << Robserved_ << '\n';
                }
            }
            else
            {
                const fileName fp = base/"Rint_vs_time.dat";
                const bool existed = isFile(fp);
                std::ofstream osf(fp.c_str(), std::ios::out | std::ios::app);
                if (osf.good())
                {
                    if (!existed) osf << "# time[s]\tRint\tRobserved\n";
                    osf << runTime.value() << '\t' << Rint_ << '\t' << Robserved_ << '\n';
                }
            }

            // Summary tables (_summary), overwritten each writeTime
            static label lastSummaryTimeIndexGlobal = -1;
            const bool newWriteTime =
                (runTime.timeIndex() != lastSummaryTimeIndexGlobal);

            if (model_ == Model::solutePermeability)
            {
                const fileName latest = sumBase/"B_sigma_latest.dat";
                if (newWriteTime)
                {
                    std::ofstream osLatest(latest.c_str(), std::ios::out | std::ios::trunc);
                    if (osLatest.good())
                        osLatest << "# time[s]\tpatch\tB\tsigmaUsed\tRobserved\n";
                    lastSummaryTimeIndexGlobal = runTime.timeIndex();
                }
                std::ofstream osLatest(latest.c_str(), std::ios::out | std::ios::app);
                if (osLatest.good())
                {
                    const scalar sigmaUsed =
                        (hasSigma_ ? sigma_ : sigmaFromUOrDefault1_(patch().index()));
                    osLatest << runTime.value() << '\t' << patch().name() << '\t'
                             << B_ << '\t' << sigmaUsed << '\t' << Robserved_ << '\n';
                }

                const fileName snap = sumBase/("B_sigma_" + runTime.timeName() + ".dat");
                const bool existedSnap = isFile(snap);
                std::ofstream osSnap(snap.c_str(), std::ios::out | std::ios::app);
                if (osSnap.good())
                {
                    if (!existedSnap)
                        osSnap << "# time[s]\tpatch\tB\tsigmaUsed\tRobserved\n";
                    const scalar sigmaUsed =
                        (hasSigma_ ? sigma_ : sigmaFromUOrDefault1_(patch().index()));
                    osSnap << runTime.value() << '\t' << patch().name() << '\t'
                           << B_ << '\t' << sigmaUsed << '\t' << Robserved_ << '\n';
                }
            }
            else
            {
                const fileName latest = sumBase/"Rint_latest.dat";
                if (newWriteTime)
                {
                    std::ofstream osLatest(latest.c_str(), std::ios::out | std::ios::trunc);
                    if (osLatest.good())
                        osLatest << "# time[s]\tpatch\tRint\tRobserved\n";
                    lastSummaryTimeIndexGlobal = runTime.timeIndex();
                }
                std::ofstream osLatest(latest.c_str(), std::ios::out | std::ios::app);
                if (osLatest.good())
                {
                    osLatest << runTime.value() << '\t' << patch().name() << '\t'
                             << Rint_ << '\t' << Robserved_ << '\n';
                }

                const fileName snap = sumBase/("Rint_" + runTime.timeName() + ".dat");
                const bool existedSnap = isFile(snap);
                std::ofstream osSnap(snap.c_str(), std::ios::out | std::ios::app);
                if (osSnap.good())
                {
                    if (!existedSnap)
                        osSnap << "# time[s]\tpatch\tRint\tRobserved\n";
                    osSnap << runTime.value() << '\t' << patch().name() << '\t'
                           << Rint_ << '\t' << Robserved_ << '\n';
                }
            }
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
