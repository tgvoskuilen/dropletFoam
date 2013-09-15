{
    Info<< "Building TEqn" << endl;
    
    Pair<tmp<volScalarField> > TSuSp = mixture.TSuSp("Cv"); //<- includes rCv, linearized evaporation enthalpy sink

    //tmp<volScalarField> rCv = mixture.rCv();
    
    /*tmp<surfaceScalarField> rhoD = fvc::interpolate(mixture.alphaVapor()) * mixture.vapor().D() + fvc::interpolate(mixture.alphaLiquid()) * mixture.liquid().D();
    
    hs = thermo.hs();
    tmp<volScalarField> test1 = fvc::laplacian(rhoD, hs)/mixture.Cp();
    tmp<volScalarField> test2 = fvc::laplacian(mixture.kByCp(turbulence->alphat()), T);
    
    Info<<"max test1 = " << Foam::max(test1()) << endl;
    Info<<"max test2 = " << Foam::max(test2()) << endl;*/
    
    tmp<volScalarField> cdT = fvc::div(mixture.DgradY(), T) / rho * runTime.deltaT();
    Info<<"Max T rise from diffusion = " << Foam::max(Foam::mag(cdT)).value() << endl;

    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T)
      + fvc::div(mixture.rhoPhi(), T)
      - fvm::Sp((fvc::ddt(rho) + fvc::div(mixture.rhoPhi())), T)
      - fvm::laplacian(mixture.kByCv(turbulence->alphat()), T)
      - fvm::div(mixture.DgradY(), T)
     ==
        mixture.rCv()*
        (
          /*  dpdt //should be DDt(p)?
          - fvc::ddt(rho, K)
          - fvc::div(mixture.rhoPhi(), K)
          + (fvc::ddt(rho) + fvc::div(mixture.rhoPhi()))*K*/
          - p*fvc::div(phi)
          + combustion->Sh()
        )
      + TSuSp.first()
      - fvm::SuSp(TSuSp.second(), T)
    );
    
    TEqn.relax();
    TEqn.solve();
    
    Info<< "T min/max   = "
        << gMin(T) << ", "
        << gMax(T) << endl;
        
        
    //Update hs with the newly calculated temperature field (Te)
    T.max(200.0);
    T.min(3500.0);
    
    mixture.setHs();
    //hs = thermo.hs();
    
    //Then update other thermo properties at new temperature
    thermo.correct();
}
