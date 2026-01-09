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


mesh.setFluxRequired(p_rgh.name());

#include "readRotationalAcceleration.H"
#include "createFvOptions.H"
#include "membraneFields.H"

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

// Effective diffusivity field Deff = DAB + nut/Sct  (auto-write if you want it dumped)
volScalarField Deff
(
    IOobject
    (
        "Deff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    DAB  // initial value; will be updated each iteration below
);

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
