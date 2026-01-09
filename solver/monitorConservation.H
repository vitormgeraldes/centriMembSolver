/*---------------------------------------------------------------------------*\
  monitorConservation.H — OpenFOAM v2506 (stationary mesh, parallel safe)
  Requisitos (fatal se faltar): CA, U, Js, Jv, CAp, phi
  Opcional: Deff (para fluxo difusivo em paredes)
  Output: <rootCase>/postProcessing/conservation/cons.dat
\*---------------------------------------------------------------------------*/

#include "fvc.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "OFstream.H"
#include "Pstream.H"
#include "OSspecific.H"

static const Foam::word DIFF_FIELD_NAME("Deff");

{
    const fvMesh& mesh_ = mesh;
    const Time&   time_ = mesh_.time();

    // ====== Require campos obrigatórios ======
    if
    (
        !mesh_.foundObject<volScalarField>("CA")  ||
        !mesh_.foundObject<volVectorField>("U")   ||
        !mesh_.foundObject<volScalarField>("Js")  ||
        !mesh_.foundObject<volScalarField>("Jv")  ||
        !mesh_.foundObject<volScalarField>("CAp") ||
        !mesh_.foundObject<surfaceScalarField>("phi")
    )
    {
        FatalErrorInFunction
            << "[monitor] Required fields missing. Need: "
            << "CA (volScalarField), U (volVectorField), "
            << "Js, Jv, CAp (volScalarField), phi (surfaceScalarField)."
            << exit(FatalError);
    }

    const volScalarField& CA  = mesh_.lookupObject<volScalarField>("CA");
    const volVectorField& U   = mesh_.lookupObject<volVectorField>("U");
    const volScalarField& JsF = mesh_.lookupObject<volScalarField>("Js");
    const volScalarField& JvF = mesh_.lookupObject<volScalarField>("Jv");
    const volScalarField& CAp = mesh_.lookupObject<volScalarField>("CAp");
    const surfaceScalarField& phi = mesh_.lookupObject<surfaceScalarField>("phi");
    (void)CAp; // apenas para forçar existência; não é usado no balanço

    const bool haveDeff = mesh_.foundObject<volScalarField>(DIFF_FIELD_NAME);
    const volScalarField* Dptr = haveDeff
        ? &mesh_.lookupObject<volScalarField>(DIFF_FIELD_NAME) : nullptr;

    // ===== Inventário (massa total) =====
    const scalarField& Vnow = mesh_.V().field();       // [m3]
    const scalarField& Cnow = CA.internalField();      // [kg/m3]
    const scalar       Mnow = gSum(Cnow * Vnow);       // [kg]

    scalar Mold = Mnow;
    const bool haveOld = (CA.timeIndex() > 0);
    if (haveOld)
    {
        const scalarField& Cold = CA.oldTime().internalField();
        Mold = gSum(Cold * Vnow);
    }

    const scalar dtStep   = max(time_.deltaTValue(), SMALL);
    const scalar dMdtStep = haveOld ? (Mnow - Mold)/dtStep : 0.0;

    static label  cache_ti   = -1;
    static scalar cache_M    = Mnow;
    static scalar cache_dMdt = dMdtStep;

    if (time_.timeIndex() != cache_ti)
    {
        cache_M    = Mnow;
        cache_dMdt = dMdtStep;
        cache_ti   = time_.timeIndex();
    }

    const scalar M    = cache_M;
    const scalar dMdt = cache_dMdt;

    // ===== Acumulação por patches =====
    scalar mIn   = 0.0;
    scalar mOut  = 0.0;
    scalar mMem  = 0.0;
    scalar mDiff = 0.0;
    scalar Qnet  = 0.0;
    scalar Qin   = 0.0;

    autoPtr<surfaceScalarField> CAfPtr
    (
        new surfaceScalarField
        (
            IOobject("CAf", time_.timeName(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE),
            fvc::interpolate(CA)
        )
    );
    const surfaceScalarField& CAf = CAfPtr();

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    forAll(pbm, pid)
    {
        const polyPatch& pp = pbm[pid];
        const word& pType = pp.type();

        // Ignorar patches acopladas (processor/cyclic/AMI) e empty/wedge
        if (pp.coupled() || pType == "empty" || pType == "wedge")
            continue;

        const scalarField& A      = mesh_.magSf().boundaryField()[pid]; // [m2]
        const scalarField& JsP    = JsF.boundaryField()[pid];           // [kg/(m2 s)]
        const scalarField& JvP    = JvF.boundaryField()[pid];           // [m/s]

        const bool memVol = (U.boundaryField()[pid].type()  == "membraneSolventFlux");
        const bool memSol = (CA.boundaryField()[pid].type() == "membraneSoluteFlux");

        // ---- Fluxo volumétrico (saída +)
        scalar Qp = 0.0;
        if (memVol)
        {
            // Usa sempre Jv (sem retrofluxo já garantido pela BC)
            Qp = gSum(JvP * A);                         // [m3/s]
        }
        else
        {
            // Convecção geral por phi (positivo = saída)
            Qp = gSum(phi.boundaryField()[pid]);        // [m3/s]
        }

        Qnet += Qp;
        if (Qp < 0) Qin += -Qp;

        // ---- Soluto pela membrana (Js) — saída +
        if (memSol)
        {
            mMem += gSum(JsP * A);                      // [kg/s]
        }

        // ---- Convecção genérica (não-membrana)
        if (!memVol && !memSol)
        {
            const scalar md =
                gSum(phi.boundaryField()[pid] * CAf.boundaryField()[pid]); // [kg/s]
            if (md >= 0) mOut += md; else mIn += -md;
        }

        // ---- Difusão: só em paredes (pType == "wall")
        if (!memVol && !memSol && haveDeff && pType == "wall")
        {
            const scalarField& Dp = Dptr->boundaryField()[pid];          // [m2/s]
            const scalarField  dC = CA.boundaryField()[pid].snGrad();    // [kg/m4]
            mDiff += - gSum(Dp * dC * A);                                // [kg/s]
        }
    }

    // Resíduos e rácios
    const scalar R    = dMdt - mIn + mOut + mMem + mDiff;          // [kg/s]
    const scalar relQ = (Qin > SMALL) ? mag(Qnet)/Qin : GREAT;

    scalar Mref = 0.0;
    Mref = max(Mref, mIn);
    Mref = max(Mref, mOut + mag(mMem) + mag(mDiff));
    Mref = max(Mref, mag(dMdt));
    const scalar relM = (Mref > SMALL) ? mag(R)/Mref : GREAT;

    const scalar CAmin = gMin(CA);
    const scalar CAmax = gMax(CA);

    if (Pstream::master())
    {
        Info<< "Conservation @ t=" << time_.timeName()
            << "  M="    << M
            << "  dM/dt="<< dMdt
            << "  mIn="  << mIn
            << "  mOut=" << mOut
            << "  mMem=" << mMem
            << "  mDiff="<< mDiff
            << "  Qnet=" << Qnet
            << "  Qin="  << Qin
            << "  relQ=" << relQ
            << "  R="    << R
            << "  relM=" << relM << nl
            << "CA[min,max] = [" << CAmin << ", " << CAmax << "]" << nl
            << endl;
    }

    // ===== Output (append) — only on solver write times, and only once (master) =====
    if (Pstream::master() && time_.outputTime())
    {
        const Time& rt = mesh_.time();

        // Case root (fora de processorN em execuções paralelas)
        const fileName outDir  = rt.globalPath()/"postProcessing"/"conservation";
        const fileName outFile = outDir/"cons.dat";
        mkDir(outDir);

        if (isDir(outDir))
        {
            OFstream os
            (
                outFile,
                IOstream::ASCII,
                IOstream::UNCOMPRESSED,
                IOstreamOption::APPEND
            );

            if (os.good())
            {
                const bool needHeader = (fileSize(outFile) == 0);
                if (needHeader)
                {
                    os<< "# time  M[kg]  dMdt[kg/s]  mIn[kg/s]  mOut[kg/s]  "
                         "mMem[kg/s]  mDiff[kg/s]  Qnet[m3/s]  Qin[m3/s]  "
                         "relQ[-]  R[kg/s]  relM[-]\n";
                }

                os<< rt.timeOutputValue() << "  "
                  << M     << "  "
                  << dMdt  << "  "
                  << mIn   << "  "
                  << mOut  << "  "
                  << mMem  << "  "
                  << mDiff << "  "
                  << Qnet  << "  "
                  << Qin   << "  "
                  << relQ  << "  "
                  << R     << "  "
                  << relM  << "\n";
            }
            else
            {
                Warning<< "[monitor] Could not open for append: " << outFile << nl;
            }
        }
        else
        {
            Warning<< "[monitor] mkDir failed? " << outDir << nl;
        }
    }
}
/*---------------------------------------------------------------------------*\
  End monitorConservation.H
\*---------------------------------------------------------------------------*/
