//File : readRotationalAcceleration.H
#include "readGravitationalAcceleration.H"   // why: provides g
#include "readhRef.H"                        // why: provides hRef
#include "readpRef.H"

Info<< "\nReading rRef" << endl;             //reference radius of the contact surface
    uniformDimensionedScalarField rRef
    (
        IOobject
        (
            "rRef",
            runTime.constant(),
            mesh.thisDb(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedScalar(word::null, dimLength, Zero)
     );


// Read SRF Properties
IOdictionary dictSRF
(
    IOobject
    (
        "SRFProperties",
        U.time().constant(),
        U.db(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

vector Rorigin(dictSRF.lookup("origin"));
vector axis (dictSRF.lookup("axis"));
scalar rpm  = readScalar(dictSRF.lookup("rpm"));

// Normalize the axis (guard against zero magnitude)
{
    const scalar m = mag(axis);
    if (m <= SMALL)
    {
        FatalErrorInFunction
            << "SRFProperties.axis has zero (or near-zero) magnitude." << nl
            << "Provide a non-zero axis vector." << exit(FatalError);
    }
    axis /= m;
}

dimensionedVector omega("omega", dimless/dimTime, vector(0,0,0));
dimensionedVector origin("origin", dimLength, vector(0,0,0));
dimensionedScalar conver("conver", dimless/dimTime, 1);
dimensionedScalar converLenght("converLenght", dimLength, 1);

// Use precise numeric π without extra includes; no API changes
conver.value() = (2.0*3.14159265358979323846/60.0)*rpm;
origin = Rorigin*converLenght;
omega  = conver*axis;

volScalarField gh
(
    IOobject
    (
        "gh",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("gh", dimensionSet(0,2,-2,0,0), 0)
);


//set ghRef - gravita Centrifugal reference potential
dimensionedScalar ghRef = (- mag(g)*hRef)
                        + 0.5*sqr(mag(omega))*sqr(rRef);

Info << "omega   = " << omega   << nl
     << "origin  = " << origin  << nl
     << "rpm     = " << rpm     << nl << endl;


Info<< "Calculating potential gravitational Centrifugal  energy field g.h\n" << endl;

// cells
gh =
    (g & mesh.C())
  + 0.5*magSqr( omega ^ (mesh.C()  - origin) )
  - ghRef;

// faces
surfaceScalarField ghf =
    (g & mesh.Cf())
  + 0.5*magSqr( omega ^ (mesh.Cf() - origin) )
  - ghRef;
