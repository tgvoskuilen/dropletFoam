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
    Info<< "min,max Uvest = " << min(mag(Uv)).value() << ", " << max(mag(Uv)).value() << endl;
    // Set phiv = phivEst
    phiv = 
    (
        (fvc::interpolate(Uv*rhov) & mesh.Sf())*alphafv
      //+ fvc::ddtPhiCorr(rAUv, rhov, Uv, phiv) //TODO: Fix near interface
    );
    
    
// -- LIQUID SETUP -------
    Info<< "Setting up liquid fields for pEqn" << endl;
    
    // Create fields    
    volScalarField rAUL(1.0/UEqnL.A());
    volScalarField rhorAUL(rhoL * rAUL);
    surfaceScalarField rhorAULf(fvc::interpolate(rhoL*rAUL));

    // Set UL = UL_est
    UL = rAUL*UEqnL.H();
    Info<< "min,max ULest = " << min(mag(UL)).value() << ", " << max(mag(UL)).value() << endl;
    // Set phiL = phiLEst
    phiL = 
    (
        (fvc::interpolate(UL*rhoL) & mesh.Sf())*alphafL
      //+ fvc::ddtPhiCorr(rAUL, rhoL, UL, phiL) //TODO: Fix near interface
    );
    
    /*phiL += 
    (
        (mixture.interface().stf())*rhorAULf*mesh.magSf()
    );*/
    
    
// -- SOLUTION LOOP ------
    Info<< "Starting pressure solution loop" << endl;
    
    tmp<volScalarField> vSu = vapor.p_Su(interface, pL, rAUv);
    tmp<volScalarField> vSp = vapor.p_Sp(interface, rAUv);
    tmp<volScalarField> LSu = liquid.p_Su(interface, pv, rhorAUL);
    tmp<volScalarField> LSp = liquid.p_Sp(interface, rhorAUL);
    
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqnv
        (
            fvm::ddt(psiv*vapor.alphaCorr(interface), pv)
          + fvc::div(phiv) // <- Really phivEst
          - fvm::laplacian(rhorAUvf*alphafv, pv)
          - vSu()
          + fvm::Sp( vSp(), pv )
        );

        fvScalarMatrix pEqnL
        (
            fvc::div(phiL) // <- Really phiLEst
          - fvm::laplacian(rhorAULf*alphafL, pL)
          - LSu()
          + fvm::Sp( LSp() , pL )
        );
        
        solve
        (
            pEqnv,
            mesh.solver(p.select(pimple.finalInnerIter()))
        );
        
        solve
        (
            pEqnL,
            mesh.solver(p.select(pimple.finalInnerIter()))
        );
        
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
    
    Info<< "min,max Uv = " << min(mag(Uv)).value() << ", " << max(mag(Uv)).value() << endl;
    Info<< "min,max UL = " << min(mag(UL)).value() << ", " << max(mag(UL)).value() << endl;
    
    Info<< "min,max pv = " << min((pv)).value() << ", " << max((pv)).value() << endl;
    Info<< "min,max pL = " << min((pL)).value() << ", " << max((pL)).value() << endl;
}
