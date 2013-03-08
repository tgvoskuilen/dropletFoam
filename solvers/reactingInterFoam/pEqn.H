{
    surfaceScalarField alphaLf(fvc::interpolate(alphaL));
    surfaceScalarField alphaVf(fvc::interpolate(alphaV));
    
    surfaceScalarField rhoLf(fvc::interpolate(rhoL));
    surfaceScalarField rhoVf(fvc::interpolate(rhoV));
      
    volScalarField rAUL(1.0/UEqnL.A());
    volScalarField rAUV(1.0/UEqnV.A());
    
    surfaceScalarField rAlphaAULf(fvc::interpolate(alphaL*rAUL));
    surfaceScalarField rAlphaAUVf(fvc::interpolate(alphaV*rAUV));
    
    // Estimate U with no force terms
    UL = rAUL*UEqnL.H();
    UV = rAUV*UEqnV.H();
    
    // Create flux estimates
    surfaceScalarField phiEstL
    (
        "phiEstL",
        (fvc::interpolate(UL) & mesh.Sf())
      + fvc::ddtPhiCorr(rAUL, rhoL, UL, phiL)
    );
    
    surfaceScalarField phiEstV
    (
        "phiEstV",
        (fvc::interpolate(UV) & mesh.Sf())
      + fvc::ddtPhiCorr(rAUV, rhoV, UV, phiV)
    );

    // Add surface tension and gravity forces to flux estimates
    phiEstL +=
    (
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rhoL)
        )*rAlphaAULf*mesh.magSf()
    );
    
    phiEstV +=
    (
        (
          - ghf*fvc::snGrad(rhoV)
        )*rAlphaAUVf*mesh.magSf()
    );
    
    //Get total estimated volume flux
    surfaceScalarField phiEst(phiEstL*alphaLf + phiEstV*alphaVf);
    
    surfaceScalarField Dp
    (
        "Dp",
        mag
        (
            alphaLf*rAlphaAULf/rhoLf
          + alphaVf*rAlphaAUVf/rhoVf
        )
    );
    
    // Dp/Dt term
    //   compressibleInterFoam uses DrhoDt here, fully explicit, rather than
    //   separate out the compressibility factor and make time-implicit.
    //   Transonic cases make this fully implicit.
    fvScalarMatrix p_DpDt = 
    (
        fvm::ddt(p)
      + fvc::div(phiV, p)
      - fvc::Sp(fvc::div(phiV), p)
    );
    
    // Solve for pressure to enforce continuity (based on phi)
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_IC
        (
            fvc::div(phiEst)
          - fvm::laplacian(Dp, p)
          - mixture.S_evap()
        );

        solve
        (
            alphaV/p*p_DpDt + p_IC,
            mesh.solver(p.select(pimple.finalInnerIter()))
        );

        if (pimple.finalNonOrthogonalIter())
        {
            surfaceScalarField mSfGradp(p_IC.flux()/Dp);
            
            phiL = phiEstL + rAlphaAULf/rhoLf*mSfGradp;
            phiV = phiEstV + rAlphaAUVf/rhoVf*mSfGradp;
            phi = alphaLf*phiL + alphaVf*phiV;

            p.relax();
            mSfGradp = p_IC.flux()/Dp;
            
            UL += fvc::reconstruct((g & mesh.Sf()) + rAlphaAULf/rhoLf*mSfGradp);
            UL.correctBoundaryConditions();
            
            UV += fvc::reconstruct((g & mesh.Sf()) + rAlphaAUVf/rhoVf*mSfGradp);
            UV.correctBoundaryConditions();
            
            U = alphaL*UL + alphaV*UV;
        }
    }
        
    //Update stored properties with new pressure
    thermo.correct();
            
    Info<< "min,max p = " << gMin(p)/1e5 << ", " 
        << gMax(p)/1e5 << " bar" << endl;
        
    Info<< "min,max U = " << min(mag(U)).value() 
        << ", " << max(mag(U)).value() << endl;
}
