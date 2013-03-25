{
    // For this to work I need to pull T out of thermo, since thermo only
    //  allows const access to T
    //  Te is the external temperature field (external to thermo)
    
    tmp<volScalarField> Sh = mixture.Sh();
    
    Info<< "Min,max Sh = " << Foam::min(Sh()).value() 
        << ", " << Foam::max(Sh()).value() << endl;
    
    //Te = thermo.T();
    
    divU = fvc::div(phi);

    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T)
      + fvc::div(mixture.rhoPhi(), T)
      - fvm::Sp(fvc::ddt(rho) + fvc::div(mixture.rhoPhi()), T)
      - fvm::laplacian(mixture.kByCv(turbulence->muEff()), T)
      + p*divU*mixture.rCv()
     ==
        Sh()*mixture.rCv()
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
