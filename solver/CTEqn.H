{
    // Deff = DAB + nut/Sct (units m^2/s)
    DTeff = DT + turbulence->nut()*(1.0/Sct);


Info << "Solve the solute continuity equation" << endl;
    fvScalarMatrix CTEqn
    (
        fvm::ddt(CT)
      + fvm::div(phi, CT)
      - fvm::laplacian(DTeff, CT)
     ==

        fvOptions(CT)
    );

    CTEqn.relax();

    CTEqn.solve();

    CT.correctBoundaryConditions();


}
