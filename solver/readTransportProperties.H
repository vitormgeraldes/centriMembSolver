// file: readTransportProperties.H

// Read your polynomial coefficients (as you already do)
IOdictionary transportPropertiesDict
(
    IOobject
    (
        "transportProperties",
        U.time().constant(),
        U.db(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Info << "Begin readTransportProperties" << endl;

// diffusivity DAB(CA) = DAB0 + ...
dimensionedScalar DAB0("DAB0", dimensionSet(0,2,-1,0,0,0,0), transportPropertiesDict.subDict("diffusivity"));
dimensionedScalar DAB1("DAB1", dimensionSet(-1,5,-1,0,0,0,0), transportPropertiesDict.subDict("diffusivity"));
dimensionedScalar DAB2("DAB2", dimensionSet(-2,8,-1,0,0,0,0), transportPropertiesDict.subDict("diffusivity"));
dimensionedScalar DAB3("DAB3", dimensionSet(-3,11,-1,0,0,0,0), transportPropertiesDict.subDict("diffusivity"));
dimensionedScalar DAB4("DAB4", dimensionSet(-4,14,-1,0,0,0,0), transportPropertiesDict.subDict("diffusivity"));

// density (if you need it elsewhere)
dimensionedScalar beta("beta", dimLength*dimLength*dimLength/dimMass, transportPropertiesDict.subDict("density"));
dimensionedScalar rho0("rho0", dimDensity, transportPropertiesDict.subDict("density"));

// DAB field (solver owns this)
volScalarField DAB
(
    IOobject
    (
        "DAB",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
    ),
    mesh,
    DAB0
);

