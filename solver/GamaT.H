// Gama.H — export Γ from membraneSoluteFlux BC into volScalarField Gama

Info<< "Compute Gama (export to volScalarField)\n";
{
    const fvMesh& mesh_ = GamaT.mesh();

    // One-time: cache pointer to CA field and membrane patch ids
    static bool init = false;
    static const volScalarField* CTfieldPtr_ = nullptr;  // pointer to field "CA"
    static List<label> memPatchIds;

    if (!init)
    {
        CTfieldPtr_ = &mesh_.lookupObject<volScalarField>("CT");

        const fvBoundaryMesh& bnd = mesh_.boundary();
        DynamicList<label> ids; ids.reserve(bnd.size());
        forAll(bnd, patchI)
        {
            const auto& CTpatch = CTfieldPtr_->boundaryField()[patchI];
            if (isA<Foam::membraneTracerFluxFvPatchScalarField>(CTpatch))
                ids.append(patchI);
        }
        memPatchIds.transfer(ids);
        init = true;
    }

    // Export Γ only on membrane patches (vectorized assignment)
    for (const label patchI : memPatchIds)
    {
        // src: membrane BC on CA
        const auto& CTpatch = CTfieldPtr_->boundaryField()[patchI];
        const auto& bc = refCast<const Foam::membraneTracerFluxFvPatchScalarField>(CTpatch);
        const scalarField& GTp = bc.GamaT();     // per-face Γ (already filled by BC)

        // dst: Gama field on the same patch
        fvPatchScalarField& pfG = GamaT.boundaryFieldRef()[patchI];

        if (GTp.size() == pfG.size())           // quick guard
        {
            pfG == GTp;                         // vectorized copy
        }
    }
}
