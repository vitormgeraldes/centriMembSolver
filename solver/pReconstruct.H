// Pressure reconstruction (unchanged)
Info << "Absolute Pressure reconstruction " << endl;
p = (p_rgh + rhok*gh)*rho0 + pRef;

p.correctBoundaryConditions();
