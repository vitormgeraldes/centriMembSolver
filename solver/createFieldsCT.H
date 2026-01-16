// createFieldsCT.H

const dictionary& cDict = runTime.controlDict();

Switch solveTracer(false);

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

Info<< "solveTracer = " << solveTracer << nl << endl;


// ------------------------------------------------------------
// CT field (real if ON, dummy if OFF)
// ------------------------------------------------------------
autoPtr<volScalarField> CTPtr;

if (solveTracer)
{
    Info<< "Reading field CT" << endl;

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
}
else
{
    CTPtr.reset
    (
        new volScalarField
        (
            IOobject
            (
                "CT",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0.0)
        )
    );
}

volScalarField& CT = CTPtr();

dimensionedScalar DT0
(
    "DT0",
    dimArea/dimTime,
    0.0
);

if (solveTracer)
{
    if (transportPropertiesDict.found("DT0"))
    {
        DT0 = dimensionedScalar("DT0", dimArea/dimTime, transportPropertiesDict);
    }
    else
    {
        WarningInFunction
            << "solveTracer=true but DT0 not found in constant/transportProperties. "
            << "Disabling tracer (solveTracer=false) and using DT0=0." << nl << endl;

        solveTracer = false;
        // DT0 stays zero
    }
}

volScalarField DT
(
    IOobject
    (
        "DT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        solveTracer ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    mesh,
    DT0
);


// ------------------------------------------------------------
// DTeff (exists always; only used when solveTracer = true)
// ------------------------------------------------------------
volScalarField DTeff
(
    IOobject
    (
        "DTeff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        solveTracer ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimArea/dimTime, 0.0)
);
