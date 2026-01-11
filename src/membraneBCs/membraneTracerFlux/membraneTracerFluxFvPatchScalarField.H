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
    membraneTracerFluxFvPatchScalarField
    ------------------------------------
    Boundary condition for tracer flux across a selective membrane interface.
    Couples diffusive and convective transport based on the
    Spiegler–Kedem / Kedem–Katchalsky framework.

    Implements three closure models:

    (1) intrinsicRejection
        JsT = Jv * CTp,  CTp = (1 - RTint) * CTw
        Film model:   k (CTw - Ci) = (RTint * Jv) * CTw
                     → CTw = [k / (k - RTint*Jv)] * Ci

        Example:
        membraneB
        {
            type    membraneTracerFlux;
            model   intrinsicRejection;
            RTint    0.9;           // trace intrinsic rejection [-]
            CTb     5.0;           // bulk reference [kg/m3] for diagnostics
            value   uniform 0;
        }

    (2) tracerPermeability  (Spiegler–Kedem closure)
        Fundamental SK equation:
            JsT = (1 - σT) * Jv * CTw + BT * (CTw - CTp)
          →  CTp = [(BT + (1 - σT)Jv) / (BT + Jv)] * CTw = αT * CTw
              αT = (BT + (1 - σT)Jv) / (BT + Jv)
            Js = Jv * CTp

        Example:
        membraneB
        {
            type    membraneTracerFlux;
            model   solutePermeability;
            BT       1e-3;          // tracer permeability [m/s]
            sigmaT   0.9;           // reflection coefficient [-] (optional)
            CTb     5.0;           // reference concentration [kg/m3]
            value   uniform 0;
        }

        Sigma selection rules:
          - if 'sigmaT' is specified → use it
          - else if the velocity BC is membraneSolventFlux → inherit σ
          - else σT = 1.0 (default: perfectly selective)

    (3) observedRejection  (inverse / target Robs)
        Infers RTint such that:
            RTobs_model = 1 - <Jv * CTp> / <Jv> / CTb
        using a safeguarded secant–bisection solver (MPI-safe).

        Example:
        membraneB
        {
            type    membraneTracerFlux;
            model   observedRejection;
            RTobs    0.80;          // target observed rejection [-]
            CTb     5.0;           // reference (feed/bulk) concentration [kg/m3]
            // Optional tuning parameters:
            // RTintMin  0.0;
            // RTintMax  0.9999;
            // tol      1e-6;
            // maxIter  50;
            value    uniform 0;
        }

    DIAGNOSTICS (automatic if CTb > 0)
    ----------------------------------
        GamaT = (CTw - CTb) / CTb                : concentration polarization [-]
        Robserved = 1 - <Jv * CTp> / <Jv> / CTb  : observed rejection [-]

    OUTPUT
    ------
        Writes Js, Jv, CTp, GamaT, and Robserved to the patch dictionary,
        and logs time series to:
            postProcessing/membraneTracerFlux/<patch>/RTint_vs_time.dat
        or (for SK model):
            postProcessing/membraneTracerFlux/<patch>/B_sigmaT_RTobserved_vs_time.dat

    UNITS
    -----
        CT, CTb [kg/m3]
        Jv [m/s]
        JsT [kg/m2/s]
        DTeff [m2/s]
        BT [m/s]
        sigmaT, RTint, RTobs, GamaT [-]

    REQUIREMENTS
    -------------
        The following volScalarFields **must exist** in the case database
        before runtime, with matching boundary sizes and correct dimensions:
            • "CTp" : permeate concentration [kg/m3]
            • "JsT"  : tracer flux [kg/m2/s]
            • "Jv"  : solvent (volume) flux [m/s]
        If any of these fields are missing or inconsistent, the solver aborts
        with a FatalError during updateCoeffs().

    NOTES
    -----
        • Parallel- and periodic-safe (uses returnReduce, master logging guard)
        • Numerically stable for small (k - RTint*Jv) via clamped denominators
        • Compatible with stationary meshes (no mesh motion)
        • Mandatory export of computed patch diagnostics (CTp, Js, Jv)
          ensures consistent post-processing and data coupling

\*---------------------------------------------------------------------------*/

#ifndef membraneTracerFluxFvPatchScalarField_H
#define membraneTracerFluxFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "volFields.H"
#include "vectorField.H"

namespace Foam
{

class membraneTracerFluxFvPatchScalarField : public mixedFvPatchScalarField
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
    scalar RTint_{0.0};

    // tracerPermeability
    scalar BT_{0.0};
    scalar sigmaT_{1.0};
    bool   hasSigmaT_{false};

    // observedRejection inputs/state
    scalar RTobs_{0.0};          // target observed rejection [-]
    scalar RTobsAchieved_{0.0};  // achieved (from inverse routine) [-]
    scalar CTb_{0.0};           // reference [kg/m3]

    // inverse controls (defaults OK; user may omit)
    scalar RTintMin_{0.0};
    scalar RTintMax_{0.9999};
    scalar tol_{1e-6};
    int    maxIter_{50};

    Model  model_{Model::intrinsicRejection};

    // ---- Per-face diagnostics ----
    mutable scalarField JsT_;     // [kg/m2/s]
    mutable scalarField Jv_;     // [m/s]
    mutable scalarField CTp_;    // [kg/m3]
    mutable scalarField GamaT_;   // [-]

    // ---- Patch-scalar diagnostics ----
    mutable scalar RTobserved_{0.0};

    // postProcessing guard (avoid duplicate lines)
    mutable label  lastLoggedTimeIndex_{-1};

    // ----------------- Helpers -----------------
    scalar sigmaTFromUOrDefault1_(const label patchI) const;

    // For a given RTint, compute flux-weighted permeate conc (LOCAL only here)
    scalar evaluateForRTint_
    (
        const scalar RTint,
        const scalarField& Ci,
        const scalarField& kTDel,
        const scalarField& Jv,
        const scalarField& Ap,
        scalar& CTpFluxAvg
    ) const;

    // Infer RTint from Robs & CTb via safeguarded secant (MPI-safe reductions)
    scalar inferRTintFromObserved_
    (
        const scalarField& Ci,
        const scalarField& kTDel,
        const scalarField& Jv,
        const scalarField& Ap,
        scalar& RTobsAchieved
    ) const;

public:

    //- Runtime type information
    TypeName("membraneTracerFlux");

    // ----------------- Constructors -----------------
    membraneTracerFluxFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    membraneTracerFluxFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    membraneTracerFluxFvPatchScalarField
    (
        const membraneTracerFluxFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    membraneTracerFluxFvPatchScalarField
    (
        const membraneTracerFluxFvPatchScalarField&
    );

    membraneTracerFluxFvPatchScalarField
    (
        const membraneTracerFluxFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    // ----------------- Accessors -----------------
    inline const scalarField& JsT()  const { return JsT_; }
    inline const scalarField& Jv()  const { return Jv_; }
    inline const scalarField& CTp() const { return CTp_; }
    inline const scalarField& GamaT()const { return GamaT_; }
    inline scalar RTobserved() const { return RTobserved_; }

    // ----------------- Core -----------------
    virtual void updateCoeffs();

    // ----------------- I/O -----------------
    virtual void write(Ostream& os) const;
};

} // End namespace Foam

#endif
// ************************************************************************* //
