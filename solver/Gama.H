// Gama.H — export Γ from membraneSoluteFlux BC into volScalarField Gama

Info<< "Compute Gama (export to volScalarField)\n";
{
    const fvMesh& mesh_ = Gama.mesh();

    // One-time: cache pointer to CA field and membrane patch ids
    static bool init = false;
    static const volScalarField* CAfieldPtr_ = nullptr;  // pointer to field "CA"
    static List<label> memPatchIds;

    if (!init)
    {
        CAfieldPtr_ = &mesh_.lookupObject<volScalarField>("CA");

        const fvBoundaryMesh& bnd = mesh_.boundary();
        DynamicList<label> ids; ids.reserve(bnd.size());
        forAll(bnd, patchI)
        {
            const auto& Cpatch = CAfieldPtr_->boundaryField()[patchI];
            if (isA<Foam::membraneSoluteFluxFvPatchScalarField>(Cpatch))
                ids.append(patchI);
        }
        memPatchIds.transfer(ids);
        init = true;
    }

    // Export Γ only on membrane patches (vectorized assignment)
    for (const label patchI : memPatchIds)
    {
        // src: membrane BC on CA
        const auto& Cpatch = CAfieldPtr_->boundaryField()[patchI];
        const auto& bc = refCast<const Foam::membraneSoluteFluxFvPatchScalarField>(Cpatch);
        const scalarField& Gp = bc.Gama();     // per-face Γ (already filled by BC)

        // dst: Gama field on the same patch
        fvPatchScalarField& pfG = Gama.boundaryFieldRef()[patchI];

        if (Gp.size() == pfG.size())           // quick guard
        {
            pfG == Gp;                         // vectorized copy
        }
    }
}
