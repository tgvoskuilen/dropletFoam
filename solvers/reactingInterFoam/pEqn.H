{
    rho = thermo.rho();
    
    // First part of density update, only changes compressible regions
    //  since psi is psiV*alphaV and is zero in the liquid regions
    //thermo.rho() -= thermo.psi() * p;
    
    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rAUf(fvc::interpolate(rAU));

    // Dp/Dt term
    //   compressibleInterFoam uses DrhoDt here, fully explicit, rather than
    //   separate out the compressibility factor and make time-implicit.
    //   Transonic cases make this fully implicit.
    fvScalarMatrix p_DpDt = 
    (
        fvm::ddt(p)
      + fvc::div(phi, p)
      - fvc::Sp(fvc::div(phi), p)
    );
    
    // Estimate U with no force terms
    U = rAU*UEqn.H();
    
    forAll(U,cellI)
    {
        U[cellI].z() = 0.0;
    }
    
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
        Pair<tmp<volScalarField> > pSuSp = mixture.SuSp_evap();
        
        fvScalarMatrix p_IC
        (
            fvc::div(phi)
          - fvm::laplacian(rAUf, p)
          - pSuSp.first()
          + fvm::SuSp(pSuSp.second(), p)
          //- mixture.S_evap()
        );

        solve
        (
            max(alphaV,0.0)*(p_DpDt/p - fvc::DDt(phi,T)/T) + p_IC,
            mesh.solver(p.select(pimple.finalInnerIter()))
        );
        
        //p.max(dimensionedScalar("pMin",p.dimensions(),1000.0));

        if (pimple.finalNonOrthogonalIter())
        {
            // Second part of thermodynamic density update
            //thermo.rho() += thermo.psi()*p;

            mixture.divComp() = (p_DpDt & p)/p - fvc::DDt(phi,T)/T;

            phi += p_IC.flux();
            //p.relax(); //from Euler solver...
        }
    }
    
    // Continuity error check checks whether thermo.rho satisfies continuity
    // with rhoPhi field
    // However, since we aren't updated rhoPhi, this would not make sense
    // so we don't check, and we instead just update thermo.rho with psi*dp
    
    
    U += rAU*fvc::reconstruct((phi - phiEst)/rAUf);
    U.correctBoundaryConditions();    
            
    //p.max(dimensionedScalar("pMin",p.dimensions(),10000.0));
    
    thermo.correct();
    
    Info<< "min,max p = " << gMin(p)/1e5 << ", " 
        << gMax(p)/1e5 << " bar" << endl;
        
    Info<< "min,max U = " << min(mag(U)).value() 
        << ", " << max(mag(U)).value() << endl;
}
