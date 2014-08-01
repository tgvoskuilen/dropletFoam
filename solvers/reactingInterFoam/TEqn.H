{
    Info<< "Building TEqn" << endl;
    
    Pair<tmp<volScalarField> > TSuSp = mixture.TSuSp(); //<- includes rCp, linearized evaporation enthalpy sink

    tmp<volScalarField> rCv = mixture.rCv();
    
    /*tmp<surfaceScalarField> rhoD = fvc::interpolate(mixture.alphaVapor()) * mixture.vapor().D() + fvc::interpolate(mixture.alphaLiquid()) * mixture.liquid().D();
    
    hs = thermo.hs();
    tmp<volScalarField> test1 = fvc::laplacian(rhoD, hs)/mixture.Cp();
    tmp<volScalarField> test2 = fvc::laplacian(mixture.kByCp(turbulence->alphat()), T);
    
    Info<<"max test1 = " << Foam::max(test1()) << endl;
    Info<<"max test2 = " << Foam::max(test2()) << endl;*/
    
    
    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T)
      + fvm::div(rhoPhi, T)
      - fvm::Sp((fvc::ddt(rho) + fvc::div(rhoPhi)), T)
      - fvm::laplacian(mixture.kByCv(turbulence->alphat()), T)
     ==
        mixture.rCv()*
        (
            fvc::div(fvc::absolute(phi, U), p)
          //- fvc::ddt(rho, K)
          //- fvc::div(rhoPhi, K)
          //+ (fvc::ddt(rho) + fvc::div(mixture.rhoPhi()))*K
          + reaction->Sh()
        )
      + TSuSp.first()
      - fvm::Sp(TSuSp.second(), T)
    );
    
    TEqn.relax();
    TEqn.solve();
    
    Info<< "T min/max   = "
        << gMin(T) << ", "
        << gMax(T) << endl;
        
        
    //Update hs with the newly calculated temperature field (Te)
    T.max(200.0);
    T.min(3500.0);
    
    mixture.setHE();
    //hs = thermo.hs();
    
    //Then update other thermo properties at new temperature
    thermo.correct();
}
