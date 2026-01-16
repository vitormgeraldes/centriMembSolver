Info<< "Reading membrane related fields for Trancer solute \n" << endl;

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

