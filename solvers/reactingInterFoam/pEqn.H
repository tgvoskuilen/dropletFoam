{
    const volScalarField& alphaV = mixture.alphaVapor();
    rho = thermo.rho();
    
    // First part of density update, only changes compressible regions
    //  since psi is psiV*alphaV and is zero in the liquid regions
    thermo.rho() -= thermo.psi() * p;
    
    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rAUf(fvc::interpolate(rAU));

    // Dp/Dt term
    //   compressibleInterFoam uses DrhoDt here, fully explicit, rather than
    //   separate out the compressibility factor.
    fvScalarMatrix p_DpDt = 
    (
        fvm::ddt(p)
      + fvc::div(phi, p)
      - fvc::Sp(fvc::div(phi), p)
    );

    // Estimate U with no force terms
    U = rAU*UEqn.H();
    
    surfaceScalarField phiEst
    (
        "phiEst",
        (fvc::interpolate(U) & mesh.Sf())
      + fvc::ddtPhiCorr(rAU, rho, U, phi)
    );

    // Add surface tension and gravity forces to flux estimate
    phi = phiEst + 
    (
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );
    
    // Solve for pressure to enforce continuity (based on phi)
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_IC
        (
            fvc::div(phi)
          - fvm::laplacian(rAUf, p)
          - mixture.S_evap()
        );

        solve
        (
            alphaV/p*p_DpDt + p_IC,
            mesh.solver(p.select(pimple.finalInnerIter()))
        );
        
        //p.max(dimensionedScalar("pMin",p.dimensions(),10000.0));

        if (pimple.finalNonOrthogonalIter())
        {
            // Second part of thermodynamic density update
            thermo.rho() += thermo.psi()*p;

            phi += p_IC.flux();
        }
    }
    
    U += rAU*fvc::reconstruct((phi - phiEst)/rAUf);
    U.correctBoundaryConditions();
    
    // Continuity error check checks whether thermo.rho satisfies continuity
    // with rhoPhi field
    // However, since we aren't updated rhoPhi, this would not make sense
    // so we don't check, and we instead just update thermo.rho with psi*dp
            
            
    p.max(dimensionedScalar("pMin",p.dimensions(),10000.0));

    
    Info<< "min,max p = " << gMin(p)/1e5 << ", " 
        << gMax(p)/1e5 << " bar" << endl;
    Info<< "min,max U = " << min(mag(U)).value() 
        << ", " << max(mag(U)).value() << endl;
}
