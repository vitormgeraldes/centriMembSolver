// Rosby.H -- Compute Rossby number field Ro
// Ro = |U| / (2 |omega| L)
// L = |C - origin| (radial distance to SRF origin)

{
    const scalar omegaMag = mag(omega).value();

    volScalarField Rosby
    (
        IOobject
        (
            "Rosby",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE   // não dependemos do runTime.write()
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    scalarField& RoI = Rosby.ref();
    const vectorField& Cc = mesh.C();

    if (omegaMag <= SMALL)
    {
        RoI = scalar(0.0);
    }
    else
    {
        forAll(RoI, i)
        {
            const scalar Uc = mag(U[i]);
            const scalar L  = mag(Cc[i] - origin.value());
            const scalar den = max(2.0*omegaMag*L, SMALL);
            RoI[i] = Uc/den;
        }
    }

    forAll(Rosby.boundaryFieldRef(), patchI)
    {
        fvPatchScalarField& RoP = Rosby.boundaryFieldRef()[patchI];
        const fvPatchVectorField& UP = U.boundaryField()[patchI];
        const vectorField& Cf = mesh.Cf().boundaryField()[patchI];

        if (omegaMag <= SMALL)
        {
            RoP == scalar(0.0);
        }
        else
        {
            forAll(RoP, faceI)
            {
                const scalar Uf = mag(UP[faceI]);
                const scalar L  = mag(Cf[faceI] - origin.value());
                const scalar den = max(2.0*omegaMag*L, SMALL);
                RoP[faceI] = Uf/den;
            }
        }
    }

    Rosby.correctBoundaryConditions();

    // >>> THIS is the key:
    if (runTime.outputTime())
    {
        Info<< "Writing Rosby field at t=" << runTime.timeName() << nl << endl;
        Rosby.write();
    }
}
