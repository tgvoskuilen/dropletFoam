{
    const volScalarField& alphaV = mixture.alphaVapor();
    
    // First part of density update, only changes compressible regions
    //thermo.rho() -= mixture.psiV() * alphaV * p;
    
    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rAUf(fvc::interpolate(rAU));

    // Dp/Dt term
    fvScalarMatrix p_DpDt = 
    (
        fvm::ddt(p)
      + fvc::div(phi, p)
      - fvc::Sp(fvc::div(phi), p)
    );

    // Estimate U
    U = rAU*UEqn.H();
    
    surfaceScalarField phiEst
    (
        "phiEst",
        (fvc::interpolate(U) & mesh.Sf())
      + fvc::ddtPhiCorr(rAU, rho, U, phi)
    );

    phi = phiEst + 
    (
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );
    
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
            //thermo.rho() += mixture.psiV()*alphaV*p;

            phi += p_IC.flux();
        }
    }
    
    U += rAU*fvc::reconstruct((phi - phiEst)/rAUf);
    U.correctBoundaryConditions();
            
            
    p.max(dimensionedScalar("pMin",p.dimensions(),10000.0));

    
    Info<< "min,max p = " << gMin(p)/1e5 << ", " << gMax(p)/1e5 << " bar" << endl;
    Info<< "min,max U = " << min(mag(U)).value() << ", " << max(mag(U)).value() << endl;
}
