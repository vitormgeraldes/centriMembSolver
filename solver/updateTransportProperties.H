// updateTransportProperties.H  (call each outer loop, BEFORE assembling UEqn)
#include "singlePhaseTransportModel.H"

Info<< "Updating transport properties (nu(CA), DAB(CA), rho)\n" << endl;

// 1) Update the laminar ν from CA via your CAViscosity model
//    (this recomputes transport->nu() internally)
laminarTransport.correct();

// 3) Update molecular diffusivity DAB = DAB(CA)
//    Updates BOTH internalField AND boundaryField from CA at the same point.
//    The boundary update is critical on the membrane patch where Cw is large
//    (concentration polarization) and DAB(Cw) > DAB(C_bulk) by a sizeable factor:
//    using DAB0 on the membrane underestimates back-diffusion and inflates Γ.
{
    auto polyDab = [&](scalar Ci) -> scalar
    {
        return max(
            DAB0.value()
          + DAB1.value()*Ci
          + DAB2.value()*Ci*Ci
          + DAB3.value()*Ci*Ci*Ci
          + DAB4.value()*Ci*Ci*Ci*Ci,
            scalar(1e-12)
        );
    };

    // internalField
    scalarField&       d = DAB.ref();
    const scalarField& c = CA.internalField();
    forAll(d, i) { d[i] = polyDab(c[i]); }

    // boundaryField — recompute from CA face values on every patch
    volScalarField::Boundary&        dBf = DAB.boundaryFieldRef();
    const volScalarField::Boundary&  cBf = CA.boundaryField();
    forAll(dBf, patchi)
    {
        scalarField&       dp = dBf[patchi];
        const scalarField& cp = cBf[patchi];
        forAll(dp, fi) { dp[fi] = polyDab(cp[fi]); }
    }
}

// 4) Update density if you’re using a Boussinesq-like relation
rhok = 1.0 + beta*CA;
rhok.correctBoundaryConditions();
Info<< "End Updating transport properties (nu(CA), DAB(CA), rho)\n" << endl;
