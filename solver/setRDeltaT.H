// setRDeltaT.H
// Compute reciprocal local time-step for LTS (localEuler) pseudo-transient
// steady-state convergence.
//
// Two contributions are combined per cell:
//   (1) Convective:  rDeltaT_conv = sum(|phi_f|) / (2 maxCo V)
//   (2) Diffusive:   rDeltaT_diff = sum(Deff_f |Sf| deltaCoeffs_f) / (maxDiNum V)
//
// rDeltaT = max(rDeltaT_conv, rDeltaT_diff, 1/maxDeltaT), then optionally
// smoothed and damped.

{
    volScalarField& rDeltaT = trDeltaT.ref();

    const dictionary& pimpleDict = pimple.dict();

    scalar maxCo
    (
        pimpleDict.getOrDefault<scalar>("maxCo", 0.8)
    );

    scalar maxDiNum
    (
        pimpleDict.getOrDefault<scalar>("maxDiNum", 5.0)
    );

    scalar rDeltaTSmoothingCoeff
    (
        pimpleDict.getOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.1)
    );

    scalar rDeltaTDampingCoeff
    (
        pimpleDict.getOrDefault<scalar>("rDeltaTDampingCoeff", 1.0)
    );

    scalar maxDeltaT
    (
        pimpleDict.getOrDefault<scalar>("maxDeltaT", GREAT)
    );

    // Damped previous field, evaluated before any update
    tmp<volScalarField> rDeltaT0_damped;

    if
    (
        rDeltaTDampingCoeff < 1
     && runTime.timeIndex() > runTime.startTimeIndex() + 1
    )
    {
        rDeltaT0_damped = (scalar(1) - rDeltaTDampingCoeff)*rDeltaT;
    }

    Info<< "LTS time scales min/max:" << endl;

    // -----------------------------------------------------------------
    // (1) Convective time scale: 1/Δt_conv = sum(|phi_f|) / (2 maxCo V)
    // -----------------------------------------------------------------
    rDeltaT.ref() =
    (
        fvc::surfaceSum(mag(phi))()()
       /((2.0*maxCo)*mesh.V())
    );

    // Limit largest time scale (smallest reciprocal time)
    rDeltaT.clamp_min(1.0/maxDeltaT);

    {
        auto limits = gMinMax(rDeltaT.primitiveField());
        limits.reset(1/(limits.max()+VSMALL), 1/(limits.min()+VSMALL));
        Info<< "    Convective = "
            << limits.min() << ", " << limits.max() << endl;
    }

    // -----------------------------------------------------------------
    // (2) Diffusive time scale (CA equation):
    //     1/Δt_diff = sum(Deff_f |Sf| deltaCoeffs_f) / (maxDiNum V)
    //     Critical for the species transport because Pe is high and the
    //     boundary-layer cells near the membrane have very small Δx.
    // -----------------------------------------------------------------
    {
        volScalarField::Internal rDeltaT_diff
        (
            fvc::surfaceSum
            (
                fvc::interpolate(Deff)
               *mesh.magSf()
               *mesh.surfaceInterpolation::deltaCoeffs()
            )()()
           /(maxDiNum*mesh.V())
        );

        rDeltaT.primitiveFieldRef() = max
        (
            rDeltaT.primitiveField(),
            rDeltaT_diff
        );

        auto limits = gMinMax(rDeltaT_diff.field());
        limits.reset(1/(limits.max()+VSMALL), 1/(limits.min()+VSMALL));
        Info<< "    Diffusive  = "
            << limits.min() << ", " << limits.max() << endl;
    }

    // Update boundary values
    rDeltaT.correctBoundaryConditions();

    // -----------------------------------------------------------------
    // Spatial smoothing — softens cell-to-cell Δt jumps that hurt
    // convergence on graded meshes (refineWallLayer near membrane).
    // -----------------------------------------------------------------
    if (rDeltaTSmoothingCoeff < 1)
    {
        fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
    }

    // -----------------------------------------------------------------
    // Damping: limit rate of change of the time scale between iterations
    // -----------------------------------------------------------------
    if (rDeltaT0_damped)
    {
        rDeltaT.clamp_min(rDeltaT0_damped());
    }

    rDeltaT.correctBoundaryConditions();

    {
        auto limits = gMinMax(rDeltaT.field());
        limits.reset(1/(limits.max()+VSMALL), 1/(limits.min()+VSMALL));
        Info<< "    Overall    = "
            << limits.min() << ", " << limits.max() << endl;
    }
}

// ************************************************************************* //
