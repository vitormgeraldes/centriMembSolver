// updateTransportProperties.H  (call each outer loop, BEFORE assembling UEqn)
#include "singlePhaseTransportModel.H"

Info<< "Updating transport properties (nu(CA), DAB(CA), rho)\n" << endl;

// 1) Update the laminar ν from CA via your CAViscosity model
//    (this recomputes transport->nu() internally)
laminarTransport.correct();

// 3) Update molecular diffusivity DAB = DAB(CA)
{
    scalarField&       d = DAB.ref();
    const scalarField& c = CA.internalField();

    forAll(d, i)
    {
        const scalar Ci = c[i];
        d[i] = max(
            DAB0.value()
          + DAB1.value()*Ci
          + DAB2.value()*Ci*Ci
          + DAB3.value()*Ci*Ci*Ci
          + DAB4.value()*Ci*Ci*Ci*Ci,
          scalar(1e-12)
        );
    }
    DAB.correctBoundaryConditions();
}

// 4) Update density if you’re using a Boussinesq-like relation
rhok = 1.0 + beta*CA;
rhok.correctBoundaryConditions();
Info<< "End Updating transport properties (nu(CA), DAB(CA), rho)\n" << endl;
