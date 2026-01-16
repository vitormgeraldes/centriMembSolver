//createFields.H

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info<< "Creating field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "p",dimensionSet(1,-1,-2,0,0,0,0),1e5)
);


Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field CA\n" << endl;
volScalarField CA
(
    IOobject
    (
        "CA",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

#include "readTransportProperties.H"

// Kinematic density for buoyancy force
volScalarField rhok
(
    IOobject
    (
        "rhok",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, 1.0)
);

// Effective diffusivity field Deff = DAB + nut/Sct  (auto-write if you want it dumped)
volScalarField Deff
(
    IOobject
    (
        "Deff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    DAB  // initial value; will be updated each iteration below
);

mesh.setFluxRequired(p_rgh.name());

#include "readRotationalAcceleration.H"
#include "createFvOptions.H"
// Remove membraneFields include - we create membrane fields directly below

Info<< "Creating turbulence model\n" << endl;

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

// Turbulent Schmidt number (dimensionless)
scalar Sct = 0.7;
if (transportPropertiesDict.found("massTransfer"))
{
    const dictionary& mt = transportPropertiesDict.subDict("massTransfer");
    Sct = mt.lookupOrDefault<scalar>("Sct", 0.7);
}
// Prevent “set but not used” when running with -postProcess
(void)Sct;


volVectorField coriolisSrc
(
    IOobject
    (
        "coriolisSrc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "zero",
        U.dimensions()/dimTime,   // [m/s] / [s] = [m/s^2]
        vector::zero
    )
);

const dictionary& controlDict = mesh.time().controlDict();

// --------------------------------------------------------------------
// Membrane-related fields (always exist)
// --------------------------------------------------------------------
Info<< "Creating membrane related fields" << endl;

// Solute flux through membrane [kg/(m^2 s)]
volScalarField Js
(
    IOobject
    (
        "Js",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Js", dimensionSet(1,-2,-1,0,0,0,0), 0.0)
);

// Volumetric flux through membrane [m/s]
volScalarField Jv
(
    IOobject
    (
        "Jv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Jv", dimensionSet(0,1,-1,0,0,0,0), 0.0)
);

// Solute concentration at the permeate side [kg/m^3]
volScalarField CAp
(
    IOobject
    (
        "CAp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("CAp", dimensionSet(1,-3,0,0,0,0,0), 0.0)
);

// Solute polarization factor at membrane [dimensionless]
volScalarField Gama
(
    IOobject
    (
        "Gama",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Gama", dimless, 0.0)
);

// --------------------------------------------------------------------
// Tracer switch and fields (created only if solveTracer=true)
// --------------------------------------------------------------------
// Read user option for solving the tracer
Switch solveTracer(false);
{
    const dictionary& cDict = runTime.controlDict();
    if (cDict.found("centriMembSolver"))
    {
        solveTracer =
            cDict.subDict("centriMembSolver")
                 .lookupOrDefault<Switch>("solveTracer", false);
    }
    else
    {
        solveTracer = cDict.lookupOrDefault<Switch>("solveTracer", false);
    }
}

Info<< "solveTracer = " << solveTracer << nl << endl;

// Pointer holders for tracer fields.  These are only allocated when
// solveTracer=true.  Otherwise they remain empty and the tracer fields
// simply do not exist.
autoPtr<volScalarField> CTPtr;
autoPtr<volScalarField> DTPtr;
autoPtr<volScalarField> DTeffPtr;
autoPtr<volScalarField> JsTPtr;
autoPtr<volScalarField> CTpPtr;
autoPtr<volScalarField> GamaTPtr;

// Default diffusivity for the tracer (will be overwritten from dictionary)
dimensionedScalar DT0("DT0", dimArea/dimTime, 0.0);

if (solveTracer)
{
    // Require DT0 in transportProperties when tracer is active
    if (!transportPropertiesDict.found("DT0"))
    {
        FatalErrorInFunction
            << "solveTracer=true but DT0 not found in constant/transportProperties." << nl
            << "Add DT0 to transportProperties or set solveTracer=false." << nl
            << exit(FatalError);
    }

    // Read DT0 from dictionary
    DT0 = dimensionedScalar("DT0", dimArea/dimTime, transportPropertiesDict);

    Info<< "Reading tracer field CT" << endl;
    CTPtr.reset
    (
        new volScalarField
        (
            IOobject
            (
                "CT",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

    // Molecular diffusivity of the tracer
    DTPtr.reset
    (
        new volScalarField
        (
            IOobject
            (
                "DT",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            DT0
        )
    );

    // Effective diffusivity of the tracer (DT + nut/SctT); updated in CTEqn.H
    DTeffPtr.reset
    (
        new volScalarField
        (
            IOobject
            (
                "DTeff",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            DT0
        )
    );

    // Membrane-related tracer outputs
    Info<< "Creating membrane related tracer fields" << endl;
    JsTPtr.reset
    (
        new volScalarField
        (
            IOobject
            (
                "JsT",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("JsT", dimensionSet(1,-2,-1,0,0,0,0), 0.0)
        )
    );

    CTpPtr.reset
    (
        new volScalarField
        (
            IOobject
            (
                "CTp",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("CTp", dimensionSet(1,-3,0,0,0,0,0), 0.0)
        )
    );

    GamaTPtr.reset
    (
        new volScalarField
        (
            IOobject
            (
                "GamaT",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("GamaT", dimless, 0.0)
        )
    );
}
