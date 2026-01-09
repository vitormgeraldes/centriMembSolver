/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  Field          | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   Operation      | Website:  www.openfoam.com                      |
|   \\  /    And            | Version:  v2506                                 |
|    \\/     Manipulation   |                                                 |
\*---------------------------------------------------------------------------*/
/*  Derived from OpenFOAM® v2506
 *  Copyright (C) 2011-2017  OpenFOAM Foundation
 *  Copyright (C) 2021       OpenCFD Ltd.
 *
 *  Modifications:
 *      by Vítor Geraldes, University of Lisbon, 2025
 *
 *  License:
 *      This file is distributed under the terms of the GNU General Public License
 *      as published by the Free Software Foundation; either version 3 of the
 *      License, or (at your option) any later version.
 *
 *  Disclaimer:
 *      Independent downstream extension compatible with OpenFOAM® v2506.
 *      Not part of, endorsed by, or affiliated with OpenCFD Ltd.,
 *      the OpenFOAM Foundation, or ESI Group.
 */
/*---------------------------------------------------------------------------*\
    membraneSoluteFluxFvPatchScalarField
    ------------------------------------
    Boundary condition for solute flux across a selective membrane interface.
    Couples diffusive and convective transport based on the
    Spiegler–Kedem / Kedem–Katchalsky framework.

    Implements three closure models:

    (1) intrinsicRejection
        Js = Jv * Cp,  Cp = (1 - Rint) * Cw
        Film model:   k (Cw - Ci) = (Rint * Jv) * Cw
                     → Cw = [k / (k - Rint*Jv)] * Ci

        Example:
        membraneB
        {
            type    membraneSoluteFlux;
            model   intrinsicRejection;
            Rint    0.9;           // intrinsic rejection [-]
            CAb     5.0;           // bulk reference [kg/m3] for diagnostics
            value   uniform 0;
        }

    (2) solutePermeability  (Spiegler–Kedem closure)
        Fundamental SK equation:
            Js = (1 - σ) * Jv * Cw + B * (Cw - Cp)
          →  Cp = [(B + (1 - σ)Jv) / (B + Jv)] * Cw = α * Cw
              α = (B + (1 - σ)Jv) / (B + Jv)
            Js = Jv * Cp

        Example:
        membraneB
        {
            type    membraneSoluteFlux;
            model   solutePermeability;
            B       1e-3;          // solute permeability [m/s]
            sigma   0.9;           // reflection coefficient [-] (optional)
            CAb     5.0;           // reference concentration [kg/m3]
            value   uniform 0;
        }

        Sigma selection rules:
          - if 'sigma' is specified → use it
          - else if the velocity BC is membraneSolventFlux → inherit σ
          - else σ = 1.0 (default: perfectly selective)

    (3) observedRejection  (inverse / target Robs)
        Infers Rint such that:
            Robs_model = 1 - <Jv * Cp> / <Jv> / CAb
        using a safeguarded secant–bisection solver (MPI-safe).

        Example:
        membraneB
        {
            type    membraneSoluteFlux;
            model   observedRejection;
            Robs    0.80;          // target observed rejection [-]
            CAb     5.0;           // reference (feed/bulk) concentration [kg/m3]
            // Optional tuning parameters:
            // RintMin  0.0;
            // RintMax  0.9999;
            // tol      1e-6;
            // maxIter  50;
            value    uniform 0;
        }

    DIAGNOSTICS (automatic if CAb > 0)
    ----------------------------------
        Γ = (Cw - CAb) / CAb                : concentration polarization [-]
        Robserved = 1 - <Jv * Cp> / <Jv> / CAb  : observed rejection [-]

    OUTPUT
    ------
        Writes Js, Jv, CAp, Gama, and Robserved to the patch dictionary,
        and logs time series to:
            postProcessing/membraneSoluteFlux/<patch>/Rint_vs_time.dat
        or (for SK model):
            postProcessing/membraneSoluteFlux/<patch>/B_sigma_Robserved_vs_time.dat

    UNITS
    -----
        CA, CAb [kg/m3]
        Jv [m/s]
        Js [kg/m2/s]
        Deff [m2/s]
        B [m/s]
        σ, Rint, Robs, Γ [-]

    REQUIREMENTS
    -------------
        The following volScalarFields **must exist** in the case database
        before runtime, with matching boundary sizes and correct dimensions:
            • "CAp" : permeate concentration [kg/m3]
            • "Js"  : solute flux [kg/m2/s]
            • "Jv"  : solvent (volume) flux [m/s]
        If any of these fields are missing or inconsistent, the solver aborts
        with a FatalError during updateCoeffs().

    NOTES
    -----
        • Parallel- and periodic-safe (uses returnReduce, master logging guard)
        • Numerically stable for small (k - Rint*Jv) via clamped denominators
        • Compatible with stationary meshes (no mesh motion)
        • Mandatory export of computed patch diagnostics (CAp, Js, Jv)
          ensures consistent post-processing and data coupling

\*---------------------------------------------------------------------------*/

#ifndef membraneSoluteFluxFvPatchScalarField_H
#define membraneSoluteFluxFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "volFields.H"
#include "vectorField.H"

namespace Foam
{

class membraneSoluteFluxFvPatchScalarField : public mixedFvPatchScalarField
{
public:

    enum class Model
    {
        intrinsicRejection,
        solutePermeability,
        observedRejection,
        Unknown
    };

private:

    // ---- User parameters (per patch) ----
    word   UName_{"U"};

    // intrinsicRejection / storage for inferred value in observedRejection
    scalar Rint_{0.0};

    // solutePermeability
    scalar B_{0.0};
    scalar sigma_{1.0};
    bool   hasSigma_{false};

    // observedRejection inputs/state
    scalar Robs_{0.0};          // target observed rejection [-]
    scalar RobsAchieved_{0.0};  // achieved (from inverse routine) [-]
    scalar CAb_{0.0};           // reference [kg/m3]

    // inverse controls (defaults OK; user may omit)
    scalar RintMin_{0.0};
    scalar RintMax_{0.9999};
    scalar tol_{1e-6};
    int    maxIter_{50};

    Model  model_{Model::intrinsicRejection};

    // ---- Per-face diagnostics ----
    mutable scalarField Js_;     // [kg/m2/s]
    mutable scalarField Jv_;     // [m/s]
    mutable scalarField CAp_;    // [kg/m3]
    mutable scalarField Gama_;   // [-]

    // ---- Patch-scalar diagnostics ----
    mutable scalar Robserved_{0.0};

    // postProcessing guard (avoid duplicate lines)
    mutable label  lastLoggedTimeIndex_{-1};

    // ----------------- Helpers -----------------
    scalar sigmaFromUOrDefault1_(const label patchI) const;

    // For a given Rint, compute flux-weighted permeate conc (LOCAL only here)
    scalar evaluateForRint_
    (
        const scalar Rint,
        const scalarField& Ci,
        const scalarField& kDel,
        const scalarField& Jv,
        const scalarField& Ap,
        scalar& CpFluxAvg
    ) const;

    // Infer Rint from Robs & CAb via safeguarded secant (MPI-safe reductions)
    scalar inferRintFromObserved_
    (
        const scalarField& Ci,
        const scalarField& kDel,
        const scalarField& Jv,
        const scalarField& Ap,
        scalar& RobsAchieved
    ) const;

public:

    //- Runtime type information
    TypeName("membraneSoluteFlux");

    // ----------------- Constructors -----------------
    membraneSoluteFluxFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    membraneSoluteFluxFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    membraneSoluteFluxFvPatchScalarField
    (
        const membraneSoluteFluxFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    membraneSoluteFluxFvPatchScalarField
    (
        const membraneSoluteFluxFvPatchScalarField&
    );

    membraneSoluteFluxFvPatchScalarField
    (
        const membraneSoluteFluxFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    // ----------------- Accessors -----------------
    inline const scalarField& Js()  const { return Js_; }
    inline const scalarField& Jv()  const { return Jv_; }
    inline const scalarField& CAp() const { return CAp_; }
    inline const scalarField& Gama()const { return Gama_; }
    inline scalar Robserved() const { return Robserved_; }

    // ----------------- Core -----------------
    virtual void updateCoeffs();

    // ----------------- I/O -----------------
    virtual void write(Ostream& os) const;
};

} // End namespace Foam

#endif
// ************************************************************************* //
