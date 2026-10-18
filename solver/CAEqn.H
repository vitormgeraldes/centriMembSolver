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
    fvOptions.constrain(CAEqn);

    CAEqn.solve();

    // Discrete solute flux of the matrix just solved (implicit convection +
    // diffusion, with the SAME phi and coefficients the equation used).
    // Used by monitorConservation.H so that the global balance is exactly
    // the one the discretised equation satisfies.
    phiCAeqn = CAEqn.flux();
    phiAtCA  = phi;      // volumetric flux the CA equation was assembled with

    fvOptions.correct(CA);
    CA.correctBoundaryConditions();


}
