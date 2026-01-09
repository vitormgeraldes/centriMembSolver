/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  Field          | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   Operation      | Website:  www.openfoam.com                      |
|   \\  /    And            | Version:  v2506                                 |
|    \\/     Manipulation   |                                                 |
\*---------------------------------------------------------------------------*/
/*  Derived from OpenFOAM® v2506
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
    membraneSolventFluxFvPatchVectorField
    -------------------------------------
    Solvent-flux membrane boundary condition (vector U) with three models.

    PHYSICS
    -------
      Spiegler–Kedem driving force:
          Jv = Ah * ( Δp - σ Δπ )
        Δp  = p_feed - p_perm
        Δπ  = π(C_feed) - π(C_perm)
      Osmotic pressure (this implementation):
          π(C) = Σ_k a_k C^k      (no temperature correction)
          with C in kg/m^3, π in Pa, and a_k in Pa·(m^3/kg)^k.

    MODELS
    ------
    (1) SpieglerKedem
        Uses user-supplied Ah and σ; π from virialCoeffs and CA/CAp.

    (2) constantFlux
        Prescribes outward-positive flux magnitude Jv [m/s].

    (3) targetAverageFlux
        Chooses an effective Ah* so that areaAverage(Jv) ≈ JvAverage
        while keeping local Jv ∝ (Δp - σΔπ). Ah* limited to [AhMin, AhMax].

    OUTPUT
    ------
        • Per-face scalar Jv_ written (outward-positive) [m/s]
        • Accessor Jv() available to solvers/post-processing
        • For SpieglerKedem, per-write snapshot:
            postProcessing/membraneSolventFlux/<patch>/permeateFlux.dat
          with columns: time[s], Jv_avg[m/s]
        • For targetAverageFlux, per-write snapshot:
            postProcessing/membraneSolventFlux/<patch>/calibratedProperties.dat

    REQUIREMENTS
    ------------
        • The following volScalarFields must exist before runtime with matching
          boundary sizes on this patch:
            - "p"   : pressure [Pa]
            - "CA"  : feed concentration [kg/m^3] (overridable via CAName)
            - "CAp" : permeate concentration [kg/m^3] (overridable via CApName)
        • Permeate-side osmotic pressure always uses π(CAp) per face.
        • Retroflux is disallowed: Jv is clamped to be non-negative.

    UNITS
    -----
        U [m/s]; Jv [m/s]; p [Pa]; π [Pa]; Ah [m/(Pa·s)];
        CA, CAp [kg/m^3]; σ [-].

        TYPICAL USAGE (boundary entries in 0/U)
        ---------------------------------------
        // 1) Spiegler–Kedem
        /*
        membrane
        {
            type            membraneSolventFlux;
            model           SpieglerKedem;

            Ah              3e-12;                 // [m/(Pa·s)]
            sigma           0.95;                  // [-]
            pPermConst      0;                     // [Pa]
            CAName          CA;                    // feed field name
            CApName         CAp;                   // permeate field name
            virialCoeffs    ( 2.5e3  1.2e-1 );     // π = a1*C + a2*C^2

            // optional moving-average reporting
            windowAverage   true;
            Twindow         10;                    // [s]

            value           uniform (0 0 0);
        }
        */

        // 2) constantFlux
        /*
        membrane
        {
            type            membraneSolventFlux;
            model           constantFlux;

            Jv              1.0e-5;                // [m/s], outward-positive
            CAName          CA;
            CApName         CAp;

            value           uniform (0 0 0);
        }
        */

        // 3) targetAverageFlux
        /*
        membrane
        {
            type            membraneSolventFlux;
            model           targetAverageFlux;

            JvAverage       8.0e-6;                // [m/s] target area-average
            AhMin           0;                     // [m/(Pa·s)]
            AhMax           1e-10;                 // [m/(Pa·s)]

            sigma           0.9;                   // [-]
            pPermConst      0;                     // [Pa]
            CAName          CA;
            CApName         CAp;
            virialCoeffs    ( 3.0e3 );             // π = a1*C

            windowAverage   true;
            Twindow         1;                    // [s]

            value           uniform (0 0 0);
        }
«

    NOTES
    -----
        • σ defaults to 1.0 if not provided.
        • virialCoeffs are interpreted as polynomial coefficients in Pa
          for C in kg/m^3 (no temperature scaling).
        • Parallel-safe reductions for area-averages; logging is master-only.
\*---------------------------------------------------------------------------*/
#ifndef membraneSolventFluxFvPatchVectorField_H
#define membraneSolventFluxFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "volFields.H"
#include "DynamicList.H"

namespace Foam
{

class membraneSolventFluxFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
public:
    enum class Model { SpieglerKedem, constantFlux, targetAverageFlux, Unknown };

private:
    // --- Cached field pointers (reset on copy/map)
    const volScalarField* pPtr_{nullptr};
    const volScalarField* CAPtr_{nullptr};
    const volScalarField* CApPtr_{nullptr};

    // --- Parameters
    scalar Ah_{0.0};        // [m/(Pa*s)]
    scalar sigma_{1.0};
    scalar JvConst_{0.0};   // [m/s]
    List<scalar> virialCoeffs_{1,1.0}; // a_k in Pa·(m^3/kg)^k

    // Permeate-side static pressure for Δp (always used)
    scalar pPermConst_{0.0};

    // targetAverageFlux mode
    scalar JvAverage_{0.0};
    scalar AhMin_{0.0};
    scalar AhMax_{VGREAT};
    mutable scalar lastAhEff_{0.0};

    // Reporting
    bool   windowAverage_{false};
    scalar Twindow_{0.0};
    mutable scalar lastJvAvgInstant_{0.0};
    mutable scalar lastJvAvgWindow_{0.0};
    mutable DynamicList<scalar> tBuf_, jvBuf_;

    mutable scalarField Jv_;  // [m/s]
    word CAName_{"CA"};
    word CApName_{"CAp"};
    Model model_{Model::SpieglerKedem};

    scalar osmoticPressure_(scalar C) const;
    void updateWindowAveraging_(const vectorField& Up) const;

public:
    TypeName("membraneSolventFlux");

    membraneSolventFluxFvPatchVectorField
    (const fvPatch&, const DimensionedField<vector,volMesh>&);
    membraneSolventFluxFvPatchVectorField
    (const fvPatch&, const DimensionedField<vector,volMesh>&, const dictionary&);
    membraneSolventFluxFvPatchVectorField
    (const membraneSolventFluxFvPatchVectorField&, const fvPatch&,
     const DimensionedField<vector,volMesh>&, const fvPatchFieldMapper&);
    membraneSolventFluxFvPatchVectorField
    (const membraneSolventFluxFvPatchVectorField&);
    membraneSolventFluxFvPatchVectorField
    (const membraneSolventFluxFvPatchVectorField&, const DimensionedField<vector,volMesh>&);

    inline scalar sigma() const { return sigma_; }
    inline const scalarField& Jv() const { return Jv_; }

    virtual void updateCoeffs() override;
    virtual void write(Ostream& os) const override;
};

} // End namespace Foam
#endif
// ************************************************************************* //
