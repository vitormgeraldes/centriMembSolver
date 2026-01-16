// GamaT.H — export Γ from membraneTracerFlux BC into volScalarField GamaT
// The tracer polarization is only available when solveTracer=true
if (solveTracer)
{
    // Access the tracer polarization field from the autoPtr
    volScalarField& GamaT = GamaTPtr();

    Info<< "Compute GamaT (export to volScalarField)\n";
    {
        const fvMesh& mesh_ = GamaT.mesh();

        // One-time: cache pointer to CT field and membrane patch ids
        static bool init = false;
        static const volScalarField* CTfieldPtr_ = nullptr;
        static List<label> memPatchIds;

        if (!init)
        {
            // lookup the tracer field (must exist because solveTracer=true)
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
            // src: membrane BC on tracer (CT)
            const auto& CTpatch = CTfieldPtr_->boundaryField()[patchI];
            const auto& bc = refCast<const Foam::membraneTracerFluxFvPatchScalarField>(CTpatch);
            const scalarField& GTp = bc.GamaT();     // per-face Γ (already filled by BC)

            // dst: GamaT field on the same patch
            fvPatchScalarField& pfG = GamaT.boundaryFieldRef()[patchI];

            if (GTp.size() == pfG.size())           // quick guard
            {
                pfG == GTp;                         // vectorized copy
            }
        }
    }
}
