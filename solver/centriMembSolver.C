/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*\
  Downstream modifications
  -------------------------
  Derived from: buoyantBoussinesqPimpleFoam (OpenFOAM v2506, ESI-OpenCFD)

  Modified by: Vítor Geraldes, University of Lisbon
  Year: 2025
  Application: centriMembSolver

  Description:
      Transient solver for centrifugal membrane processes with optional
      buoyancy, turbulence, and steady rotational effects.

  Disclaimer:
      This independent downstream work is not part of, endorsed by,
      or affiliated with OpenCFD Ltd., the OpenFOAM Foundation, or ESI Group.

  Trademarks:
      OpenFOAM® and OpenCFD® are registered trademarks of OpenCFD Ltd.
\*---------------------------------------------------------------------------*/


#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "membraneSoluteFluxFvPatchScalarField.H"
#include "membraneTracerFluxFvPatchScalarField.H"
#include "membraneSolventFluxFvPatchVectorField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for buoyant, turbulent flow"
        " of incompressible fluids, in a single uniformly rotating frame."
        " Uses the Boussinesq approximation."
        " Uses membrane boundary conditions under uniform rotation."
        " Static mesh (no mesh motion or topology changes)."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"

    #include "createMesh.H"

    #include "createFields.H"

    // PIMPLE controller (v2506)
    #include "createPimpleControl.H"

    #include "CourantNo.H"

    // Time controls: declare & read (v2506 expects these before setInitialDeltaT)
    #include "createTimeControls.H"
    #include "readTimeControls.H"
    #include "setInitialDeltaT.H"

    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    #include "updateTransportProperties.H"

    while (runTime.run())
    {

        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"
            #include "CAEqn.H"
            #include "updateTransportProperties.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }

        }
        #include "CTEqn.H"

        #include "pReconstruct.H"

        #include "monitorConservation.H"

        #include "Gama.H"
        #include "GamaT.H"

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
