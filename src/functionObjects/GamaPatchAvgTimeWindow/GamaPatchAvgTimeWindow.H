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
 *      by Vítor Geraldes, University of Lisbon, 2025-2026
 *
 *  License: GPLv3 (or later)
 */
/*---------------------------------------------------------------------------*\
    GamaPatchAvgTimeWindow
    -----------------------
    Function object for computing:
      (i)  instantaneous area-weighted patch mean, meanS(t)
      (ii) instantaneous area-weighted spatial coefficient of variation on the patch,
           CVs(t) = stdS(t)/meanS(t)
      (iii) time-windowed (moving) mean of meanS(t) and CVs(t)

    The time window is an integration window of length "timeWindow" using
    time-step weighting (dt).

    Outputs:
      - postProcessing/<outputName>/<fileName>
      - postProcessing/KPI/Gama_<patch>.txt   (KPI = winMean(CVs))

    Columns:
      Time    meanS    stdS    CVs    winMean(meanS)    winMean(CVs)
\*---------------------------------------------------------------------------*/

#ifndef GamaPatchAvgTimeWindow_H
#define GamaPatchAvgTimeWindow_H

#include "fvMeshFunctionObject.H"
#include "volFieldsFwd.H"
#include "OFstream.H"
#include <deque>

namespace Foam
{

class GamaPatchAvgTimeWindow
:
    public functionObjects::fvMeshFunctionObject
{
    struct Sample
    {
        scalar dt;
        scalar value;
    };

    // --- dictionary parameters ---
    word   patchName_;
    word   fieldName_;
    word   outputName_;
    word   fileName_;
    scalar timeWindow_;
    Switch masterOnly_;

    // --- internal state ---
    label  patchi_;

    // moving window means (time-weighted)
    scalar Wm_, Sm_;           // for meanS
    scalar Wcv_, Scv_;         // for CVs
    std::deque<Sample> bufferMean_;
    std::deque<Sample> bufferCV_;

    // --- file output ---
    autoPtr<OFstream> out_;
    fileName outPath_;
    fileName kpiPath_;

    // --- last instantaneous values (computed in execute) ---
    scalar lastInstMeanS_;
    scalar lastInstStdS_;
    scalar lastInstCVs_;

    // --- helpers ---
    void patchAreaMeanStd
    (
        const volScalarField& f,
        scalar& mean,
        scalar& stdS
    ) const;

    void pushSample(std::deque<Sample>& buf, scalar& W, scalar& S, scalar dt, scalar value);
    void trimWindow(std::deque<Sample>& buf, scalar& W, scalar& S);

public:
    TypeName("GamaPatchAvgTimeWindow");

    GamaPatchAvgTimeWindow
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    virtual bool execute();
    virtual bool write();
};

} // namespace Foam

#endif
