
    //Load field references
    const PLICInterface& interface = mixture.interface();
    
    phase& vapor = mixture.vapor();
    const volScalarField& alphav = mixture.alphaVapor();
    surfaceScalarField& phiv = vapor.phi();
    volVectorField& Uv = vapor.U();
    tmp<surfaceScalarField> talphafv = interface.alphaVaporf();
    const surfaceScalarField& alphafv = talphafv();
    const volScalarField& rhov = vapor.rho();
    volScalarField& pv = vapor.p();
    const volScalarField& muv = vapor.mu();
    
    phase& liquid = mixture.liquid();
    const volScalarField& alphaL = mixture.alphaLiquid();
    surfaceScalarField& phiL = liquid.phi();
    volVectorField& UL = liquid.U();
    const surfaceScalarField& alphafL = interface.alphaLiquidf();
    const volScalarField& rhoL = liquid.rho();
    volScalarField& pL = liquid.p();
    const volScalarField& muL = liquid.mu();
    
    
    // Make U equations
    Info<< "Building UEqns" << endl;
    
    word divScheme("div(phi,U)");
    
    fvVectorMatrix UEqnL
    (
        fvm::ddt(rhoL*liquid.alphaCorr(interface), UL)
      + fvm::div(phiL, UL, divScheme)
      - fvm::laplacian(fvc::interpolate(muL)*alphafL, UL)
      //- fvc::div(muL*dev2(T(fvc::grad(UL)))) //TODO: Fix near interface?
      - liquid.U_Su(interface, Uv, muv)
      + fvm::Sp( liquid.U_Sp(interface), UL )
      - rhoL*alphaL*g
    );

    fvVectorMatrix UEqnV
    (
        fvm::ddt(rhov*vapor.alphaCorr(interface), Uv)
      + fvm::div(phiv, Uv, divScheme)
      //+ turbulence->divDevRhoReff(Uv) //causes crash, UEqnV.H -> 1e300
      - vapor.U_Su(interface, UL, muL)
      + fvm::Sp( vapor.U_Sp(interface), Uv )
      - rhov*alphav*g
    );
    

    UEqnL.relax();
    UEqnV.relax();
    
    if (pimple.momentumPredictor())
    {
        Info<< "Doing momentum predictor step" << endl;
        
        /*Foam::solve
        (
            UEqnL
         ==
            fvc::reconstruct
            (
                (
                  - fvc::snGrad(pL)
                ) * mesh.magSf() * alphaf
                
                //- grad p from interface
            )
        );*/
        
        /*Foam::solve
        (
            UEqnV
         ==
            fvc::reconstruct
            (
                (
                  - fvc::snGrad(pV)
                ) * mesh.magSf() * alphaf
                
                //- grad p from interface
            )
        );*/
    }
    
