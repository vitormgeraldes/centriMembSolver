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
