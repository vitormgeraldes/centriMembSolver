// CTEqn.H
// The tracer equation is solved only when solveTracer is true.
if (solveTracer)
{
    // Alias convenience references from autoPtrs
    volScalarField& CT    = CTPtr();
    volScalarField& DT    = DTPtr();
    volScalarField& DTeff = DTeffPtr();

    // Update effective diffusivity of the tracer: DTeff = DT + nut/Sct
    DTeff = DT + turbulence->nut()*(1.0/Sct);

    Info<< "Solve the tracer continuity equation (CT)" << endl;

    fvScalarMatrix CTEqn
    (
        fvm::ddt(CT)
      + fvm::div(phi, CT)
      - fvm::laplacian(DTeff, CT)
     ==
        fvOptions(CT)
    );

    CTEqn.relax();
    fvOptions.constrain(CTEqn);   // recommended when using fvOptions

    CTEqn.solve();

    fvOptions.correct(CT);
    CT.correctBoundaryConditions();
}
