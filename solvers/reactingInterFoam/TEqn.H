{
    // For this to work I need to pull T out of thermo, since thermo only
    //  allows const access to T
    //  Te is the external temperature field (external to thermo)
    
    tmp<volScalarField> Sh = mixture.Sh();
    
    Info<< "Min,max Sh = " << Foam::min(Sh()).value() 
        << ", " << Foam::max(Sh()).value() << endl;
    
    Te = thermo.T();
    
    divU = fvc::div(phi);

    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, Te)
      + fvc::div(mixture.rhoPhi(), Te)
      - fvm::laplacian(mixture.kByCv(turbulence->muEff()), Te)
      + p*divU/mixture.Cv()
     ==
       Sh()/mixture.Cv()
    );

    TEqn.relax();
    TEqn.solve();
    
    Info<< "Te min/max   = "
        << gMin(Te) << ", "
        << gMax(Te) << endl;
        
        
    //Update hs with the newly calculated temperature field (Te)
    Te.max(200.0);
    Te.min(3500.0);
    
    mixture.setHs(Te);
    hs = thermo.hs();
    
    //Then update thermo.T based on new hs (so it should match Te)
    thermo.correct(); //updates thermo.rho(p,T,alphas,Ys)
    
    volScalarField dT = Foam::mag(Te - thermo.T());
    
    Info<< "dT min/max   = "
        << gMin(dT) << ", "
        << gMax(dT) << endl;
        
    Info<< "T min/max   = "
        << gMin(T) << ", "
        << gMax(T) << endl;

    
}
