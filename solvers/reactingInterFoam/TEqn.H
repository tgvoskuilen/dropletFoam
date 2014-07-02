{
    Info<< "Building TEqn" << endl;
    
    Pair<tmp<volScalarField> > TSuSp = mixture.TSuSp("Cv"); //<- includes rCv, linearized evaporation enthalpy sink

    // Build governing equation for temperature    
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
