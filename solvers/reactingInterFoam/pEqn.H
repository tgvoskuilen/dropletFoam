{

    
// -- VAPOR SETUP -------
    Info<< "Setting up vapor fields for pEqn" << endl;
    
    // Load/create fields
    tmp<volScalarField> tpsiv = vapor.psi();
    const volScalarField& psiv = tpsiv();
    
    volScalarField rAUv(1.0/UEqnV.A());
    surfaceScalarField rhorAUvf(fvc::interpolate(rhov*rAUv));

    // Set Uv = Uv_est
    Uv = rAUv*UEqnV.H();
    
    // Set phiv = phivEst
    phiv = 
    (
        (fvc::interpolate(Uv*rhov) & mesh.Sf())*alphafv
      //+ fvc::ddtPhiCorr(rAU, rho, U, phi) //TODO: Fix near interface
    );
    
    
// -- LIQUID SETUP -------
    Info<< "Setting up liquid fields for pEqn" << endl;
    
    // Create fields    
    volScalarField rAUL(1.0/UEqnL.A());
    volScalarField rhorAUL(rhoL * rAUL);
    surfaceScalarField rhorAULf(fvc::interpolate(rhoL*rAUL));

    // Set UL = UL_est
    UL = rAUL*UEqnL.H();
    
    // Set phiL = phiLEst
    phiL = 
    (
        (fvc::interpolate(UL*rhoL) & mesh.Sf())*alphafL
      //+ fvc::ddtPhiCorr(rAU, rho, U, phi) //TODO: Fix near interface
    );
    
    
// -- SOLUTION LOOP ------
    Info<< "Starting pressure solution loop" << endl;
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqnv
        (
            fvm::ddt(psiv*vapor.alphaCorr(interface), pv)
          + fvc::div(phiv) // <- Really phivEst
          - fvm::laplacian(rhorAUvf*alphafv, pv)
          - vapor.p_Su(interface, pL, rAUv)              //<- rAUv not used
          + fvm::Sp( vapor.p_Sp(interface, rAUv), pv )   //<- rAUv not used
        );
        
        fvScalarMatrix pEqnL
        (
            fvc::div(phiL) // <- Really phiLEst
          - fvm::laplacian(rhorAULf*alphafL, pL)
          - liquid.p_Su(interface, pv, rhorAUL)
          + fvm::Sp( liquid.p_Sp(interface, rhorAUL), pL )
        );

        solve
        (
            pEqnv,
            mesh.solver(p.select(pimple.finalInnerIter()))
        );
        
        Info<< "min,max pv = " << gMin(pv)/1e5 << ", " << gMax(pv)/1e5 << " bar" << endl;
        
        solve
        (
            pEqnL,
            mesh.solver(p.select(pimple.finalInnerIter()))
        );
        
        Info<< "min,max pL = " << gMin(pL)/1e5 << ", " << gMax(pL)/1e5 << " bar" << endl;

        if (pimple.finalNonOrthogonalIter())
        {
            phiv += pEqnv.flux();
            phiL += pEqnL.flux();
            
            // Correct velocity fields
            Uv += rAUv*fvc::reconstruct(pEqnv.flux()/rhorAUvf); //TODO: Fix near interface
            Uv.correctBoundaryConditions();
              
            UL += rAUL*fvc::reconstruct(pEqnL.flux()/rhorAULf); //TODO: Fix near interface
            UL.correctBoundaryConditions(); 
        }
    }            
            
    pv.max(dimensionedScalar("pMin",pv.dimensions(),10000.0));
    pL.max(dimensionedScalar("pMin",pL.dimensions(),10000.0));
    
    Info<< "min,max pv = " << gMin(pv)/1e5 << ", " << gMax(pv)/1e5 << " bar" << endl;
    Info<< "min,max Uv = " << min(mag(Uv)).value() << ", " << max(mag(Uv)).value() << endl;
    
    Info<< "min,max pL = " << gMin(pL)/1e5 << ", " << gMax(pL)/1e5 << " bar" << endl;
    Info<< "min,max UL = " << min(mag(UL)).value() << ", " << max(mag(UL)).value() << endl;
}
