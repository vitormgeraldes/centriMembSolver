Info<< "Reading membrane related fields \n" << endl;
//#include "fvcVolumeIntegrate.H"


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
volScalarField JsT
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
);
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

volScalarField CAp
(
    IOobject("CAp", runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("CAp", dimensionSet(1,-3,0,0,0,0,0), 0.0) // kg/m3
);
volScalarField CTp
(
    IOobject("CTp", runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("CTp", dimensionSet(1,-3,0,0,0,0,0), 0.0) // kg/m3
);

volScalarField Gama
(
    IOobject("Gama", runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Gama", dimensionSet(0,0,0,0,0,0,0), 0.0) // kg/m3
);
volScalarField GamaT
(
    IOobject("GamaT", runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("GamaT", dimensionSet(0,0,0,0,0,0,0), 0.0) // kg/m3
);
