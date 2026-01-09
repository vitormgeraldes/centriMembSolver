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
    GamaPatchAvgTimeWindow
    -----------------------
    Function object for computing the time-windowed, area-weighted average of a
    scalar field over a specified boundary patch.

    Maintains a moving integration window of length “timeWindow”, continuously
    updating the average and variance using time-weighted accumulation.
    Outputs are written to:
        • postProcessing/<outputName>/areaAvg_timeWindow.dat
        • postProcessing/KPI/Gama_<patch>.txt

    KEY FEATURES
    -------------
      • Parallel-safe (uses MPI reductions, no blocking barriers)
      • Master-only or per-rank output selectable
      • Automatic directory creation under postProcessing/
      • Reports instantaneous, window-mean, variance, and standard deviation
      • Simple KPI file for external monitoring or control feedback

    EXAMPLE (system/controlDict)
    ----------------------------
        functions
        {
            GamaPatchAvgTimeWindow1
            {
                type        GamaPatchAvgTimeWindow;
                libs        ("libGamaPatchAvgTimeWindow.so");

                patchName   membraneB;        // target patch
                fieldName   Gama;             // volScalarField to average
                timeWindow  15;               // window width [s]

                masterOnly  true;             // only master writes
                writeControl writeTime;
                writeInterval 1;

                outputName  Gama_membrane_avgWin;
                fileName    areaAvg_timeWindow.dat;
            }
        }

    OUTPUT FORMAT
    --------------
        # Patch: membraneB
        # Field: Gama
        # timeWindow(s): 15
        # Columns: Time   areaAverage(Gama@membraneB)   winMean   winVar   winStd

    UNITS
    -----
        Field value  : same as input volScalarField
        Time         : [s]

    NOTES
    -----
        • Uses gSum() for distributed patch integration
        • No dependency on collated fileHandler
        • Safe under decomposePar / reconstructPar workflows
        • Suitable for steady-state or transient runs

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
    scalar W_, S1_, S2_;
    std::deque<Sample> buffer_;

    // --- file output ---
    autoPtr<OFstream> out_;
    fileName outPath_;
    fileName kpiPath_;
    bool wroteHeader_;

    // --- last instantaneous global average (from execute) ---
    scalar lastInstAvg_;

    // --- internal helpers ---
    scalar patchAreaAverage(const volScalarField& f) const;
    void pushSample(scalar dt, scalar value);
    void trimWindow();

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
