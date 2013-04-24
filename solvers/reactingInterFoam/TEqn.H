{
    Info<< "Building TEqn" << endl;
        
    divU = fvc::div(phi);
    
    Pair<tmp<volScalarField> > TSuSp = mixture.TSuSp(); //<- includes rCv, linearized evaporation enthalpy sink

    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T)
      + fvc::div(mixture.rhoPhi(), T)
      - fvm::Sp(fvc::ddt(rho) + fvc::div(mixture.rhoPhi()), T)
      - fvm::laplacian(mixture.kByCv(turbulence->muEff()), T)
      + p*divU*mixture.rCv()
     ==
        combustion->Sh()*mixture.rCv()
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
    
    mixture.setHs();
    //hs = thermo.hs();
    
    //Then update other thermo properties at new temperature
    thermo.correct();
}
