// Solve the momentum equation
Info<< "Solving momentum equation (mid-point Coriolis, no ramp)" << endl;


coriolisSrc =
    -2.0*(omega ^ U);

fvVectorMatrix UEqn
(
    fvm::ddt(U)
  + fvm::div(phi, U)
  + turbulence->divDevReff(U)
  - coriolisSrc
  ==
    fvOptions(U)
);

UEqn.relax();
fvOptions.constrain(UEqn);


// =====================================================
// 4) Momentum predictor
// =====================================================

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            (
              - ghf*fvc::snGrad(rhok)
              - fvc::snGrad(p_rgh)
            )*mesh.magSf()
        )
    );

    fvOptions.correct(U);
}
