// Solve the momentum equation
Info<< "Solving momentum equation (mid-point Coriolis, no ramp)" << endl;

// =====================================================
// 1) Mid-point velocity U^{n+1/2} ≈ (U^n + U^{n+1,m})/2
// =====================================================

volVectorField Umid
(
    IOobject
    (
        "Umid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    0.5*(U.oldTime() + U)
);

// =====================================================
// 2) Coriolis acceleration based on mid-point velocity
//     f_C = -2 Ω × U_mid
// =====================================================

coriolisSrc =
    -2.0*(omega ^ Umid);

// (Opcional) diagnóstico para ver se está a atualizar mesmo
// Info<< "max|coriolisSrc| = "
//     << gMax(mag(coriolisSrc))
//     << "  time = " << runTime.timeName()
//     << "  outer = " << pimple.corr()
//     << nl << endl;


// =====================================================
// 3) Montagem da equação de momento
// =====================================================

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
