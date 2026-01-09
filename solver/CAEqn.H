{
    // Deff = DAB + nut/Sct (units m^2/s)
    Deff = DAB + turbulence->nut()*(1.0/Sct);


Info << "Solve the solute continuity equation" << endl;
    fvScalarMatrix CAEqn
    (
        fvm::ddt(CA)
      + fvm::div(phi, CA)
      - fvm::laplacian(Deff, CA)
     ==

        fvOptions(CA)
    );

    CAEqn.relax();

    CAEqn.solve();

    CA.correctBoundaryConditions();


}
