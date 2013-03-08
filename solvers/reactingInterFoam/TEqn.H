volScalarField SpL = (1.0 - liquid.cellMask())*liquid.rho(p)/runTime.deltaT();
volScalarField SpV = (1.0 - vapor.cellMask())*vapor.rho(p)/runTime.deltaT();

{
    word divScheme("div(rho*phi,T)");

    fvScalarMatrix TEqnL
    (
        fvm::ddt(rhoL, TL)
      + fvc::div(rhoPhiL, TL, divScheme)
      - fvm::Sp(fvc::div(rhoPhiL) - liquid.m_evap_sum(), TL)
      - fvm::laplacian(liquid.kByCp(), TL)
      //+ p*divU*mixture.rCv()
     ==
        SpL*TV
      - fvm::Sp(SpL, TL)
     //  Sh()*mixture.rCv()
    );
    
    fvScalarMatrix TEqnV
    (
        fvm::ddt(rhoV, TV)
      + fvc::div(rhoPhiV, TV, divScheme)
      - fvm::Sp(fvc::div(rhoPhiV) - vapor.m_evap_sum(), TV)
      - fvm::laplacian(vapor.kByCp() + turbulence->muEff(), TV)
      //+ p*divU*mixture.rCv()
     ==
        SpV*TL
      - fvm::Sp(SpV, TV)
     //  Sh()*mixture.rCv()
    );

    TEqnL.relax();
    TEqnL.solve();
    
    TEqnV.relax();
    TEqnV.solve();
    
    Info<< "TL min/max   = "
        << gMin(TL) << ", "
        << gMax(TL) << endl;
        
    Info<< "TV min/max   = "
        << gMin(TV) << ", "
        << gMax(TV) << endl;
            
    //Combine TL and TV into thermo.T and update hs accordingly
    mixture.setHs();
    
    //Then update other thermo properties at new temperature
    thermo.correct();
}
