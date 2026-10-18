/*---------------------------------------------------------------------------*\
  monitorConservation.H — global solute balance CONSISTENT with the discretised
  CA equation (v1.2.0, 18-Oct-2026).

  The balance uses the face fluxes of the CA matrix that was just solved
  (phiCAeqn = CAEqn.flux(): implicit convection + diffusion, with the phi the
  equation was assembled with) and the time derivative of the same scheme
  (fvc::ddt). The residual

      R = sum_cells V*ddt(CA) + sum_boundary phiCAeqn - mSp

  is therefore zero up to linear-solver tolerance, whatever the boundary
  conditions. mSp is the cell source of the 'bounded' convection scheme,
  -CA*div(phi), zero when continuity is satisfied.

  Open boundaries (fixed pressure + inletOutlet) can carry flow in BOTH
  directions through the same patch once buoyant plumes reach them. The solute
  leaving at the local concentration and replaced by feed at C0 is a genuine
  boundary flux of the equation, NOT a conservation error. It is reported per
  patch in postProcessing/conservation/exchange.dat:
      Qin, Qout       volumetric flow entering / leaving the patch [m3/s]
      mFeedIn         solute entering with the inflow faces      [kg/s]
      mExchOut        solute leaving through the outflow faces   [kg/s]
  so that the user can compute the exchange loss mExchOut - Qout*C0.

  Output: postProcessing/conservation/cons.dat (columns unchanged except mDiff,
  replaced by mSp) and exchange.dat.
\*---------------------------------------------------------------------------*/

#include "fvc.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "OFstream.H"
#include "Pstream.H"
#include "OSspecific.H"

{
    const fvMesh& mesh_ = mesh;
    const Time&   time_ = mesh_.time();

    // ===== Inventory =====
    const scalarField& Vnow = mesh_.V().field();
    const scalarField& Cnow = CA.internalField();
    const scalar Mnow = gSum(Cnow*Vnow);                                  // [kg]

    // Time derivative of the ddt scheme actually used (Euler, CN, backward, LTS)
    const scalar dMdt = gSum(fvc::ddt(CA)().primitiveField()*Vnow);      // [kg/s]

    // Cell source of the 'bounded' convection scheme: -Sp(div(phi)) -> +CA*div(phi) on the flux side
    const scalar mSp = gSum(Cnow*fvc::div(phiAtCA)().primitiveField()*Vnow);

    // ===== Boundary fluxes from the matrix (out positive) =====
    scalar mIn = 0, mOut = 0, mMem = 0, Qnet = 0, Qin = 0;
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    // per-patch exchange diagnostics
    DynamicList<word>   xName;
    DynamicList<scalar> xQin, xQout, xFeed, xExch;

    forAll(pbm, pid)
    {
        const polyPatch& pp = pbm[pid];
        if (pp.coupled() || pp.type() == "empty" || pp.type() == "wedge") continue;

        const scalarField& fp   = phiCAeqn.boundaryField()[pid];        // [kg/s] matrix flux
        const scalarField& php  = phiAtCA.boundaryField()[pid];         // [m3/s]
        const scalarField& CAb  = CA.boundaryField()[pid];              // face values used by the matrix

        const scalar mp = gSum(fp);
        const scalar Qp = gSum(php);
        Qnet += Qp; if (Qp < 0) Qin -= Qp;

        const bool memSol = (CA.boundaryField()[pid].type() == "membraneSoluteFlux");
        if (memSol)
        {
            mMem += mp;
        }
        else
        {
            if (mp >= 0) mOut += mp; else mIn -= mp;

            scalar qi = 0, qo = 0, mi = 0, mo = 0;
            forAll(php, i)
            {
                if (php[i] < 0) { qi -= php[i]; mi -= php[i]*CAb[i]; }
                else            { qo += php[i]; mo += php[i]*CAb[i]; }
            }
            reduce(qi, sumOp<scalar>()); reduce(qo, sumOp<scalar>());
            reduce(mi, sumOp<scalar>()); reduce(mo, sumOp<scalar>());
            xName.append(pp.name()); xQin.append(qi); xQout.append(qo); xFeed.append(mi); xExch.append(mo);
        }
    }

    const scalar R    = dMdt - mIn + mOut + mMem - mSp;                   // [kg/s]
    const scalar relQ = (Qin > SMALL) ? mag(Qnet)/Qin : GREAT;
    scalar Mref = max(mIn, max(mOut + mag(mMem), mag(dMdt)));
    const scalar relM = (Mref > SMALL) ? mag(R)/Mref : GREAT;

    const scalar CAmin = gMin(CA);
    const scalar CAmax = gMax(CA);

    if (Pstream::master())
    {
        Info<< "Conservation @ t=" << time_.timeName()
            << "  M=" << Mnow << "  dM/dt=" << dMdt
            << "  mIn=" << mIn << "  mOut=" << mOut << "  mMem=" << mMem
            << "  mSp=" << mSp << "  Qnet=" << Qnet << "  Qin=" << Qin
            << "  relQ=" << relQ << "  R=" << R << "  relM=" << relM << nl;
        forAll(xName, k)
        {
            if (xQout[k] > SMALL)
            {
                Info<< "Exchange @ patch " << xName[k]
                    << "  Qin=" << xQin[k] << "  Qout=" << xQout[k]
                    << "  recirculation Qout/Qnet=" << xQout[k]/max(xQin[k]-xQout[k], SMALL)
                    << "  mFeedIn=" << xFeed[k] << "  mExchOut=" << xExch[k] << nl;
            }
        }
        Info<< "CA[min,max] = [" << CAmin << ", " << CAmax << "]" << nl << endl;
    }

    if (Pstream::master() && time_.outputTime())
    {
        const fileName outDir = time_.globalPath()/"postProcessing"/"conservation";
        mkDir(outDir);
        {
            const fileName outFile = outDir/"cons.dat";
            OFstream os(outFile, IOstream::ASCII, IOstream::UNCOMPRESSED, IOstreamOption::APPEND);
            if (os.good())
            {
                if (fileSize(outFile) == 0)
                {
                    os<< "# time  M[kg]  dMdt[kg/s]  mIn[kg/s]  mOut[kg/s]  mMem[kg/s]  "
                         "mSp[kg/s]  Qnet[m3/s]  Qin[m3/s]  relQ[-]  R[kg/s]  relM[-]\n";
                }
                os<< time_.timeOutputValue() << "  " << Mnow << "  " << dMdt << "  "
                  << mIn << "  " << mOut << "  " << mMem << "  " << mSp << "  "
                  << Qnet << "  " << Qin << "  " << relQ << "  " << R << "  " << relM << "\n";
            }
        }
        {
            const fileName outFile = outDir/"exchange.dat";
            OFstream os(outFile, IOstream::ASCII, IOstream::UNCOMPRESSED, IOstreamOption::APPEND);
            if (os.good())
            {
                if (fileSize(outFile) == 0)
                {
                    os<< "# time  patch  Qin[m3/s]  Qout[m3/s]  mFeedIn[kg/s]  mExchOut[kg/s]\n";
                }
                forAll(xName, k)
                {
                    os<< time_.timeOutputValue() << "  " << xName[k] << "  " << xQin[k] << "  "
                      << xQout[k] << "  " << xFeed[k] << "  " << xExch[k] << "\n";
                }
            }
        }
    }
}
