volScalarField SpL = (1.0 - liquid.cellMask())*liquid.rho(p)/runTime.deltaT();
volScalarField SpV = (1.0 - vapor.cellMask())*vapor.rho(p)/runTime.deltaT();

{
    word divScheme("div(rho*phi,T)");
    word KdivScheme("div(rho*phi,K)");
    
    volScalarField heatTransferCoeff = 
    (
        alphaL*liquid.cellMask()*liquid.k()
      + alphaV*vapor.cellMask()*vapor.k()
    ) * mag(fvc::grad(alphaL));
    
    heatTransferCoeff.dimensionedInternalField() /= Foam::pow(mesh.V(),1.0/3.0);
    
    
    fvScalarMatrix TEqnL
    (
        fvm::ddt(rhoL, TL)
      + fvm::div(rhoPhiL, TL, divScheme)
      - fvm::Sp(fvc::ddt(rhoL) + fvc::div(rhoPhiL) - liquid.m_evap_sum(), TL)
      - fvm::laplacian(liquid.kByCp(), TL)
     ==
        heatTransferCoeff*TV/liquid.Cp()
      - fvm::Sp(heatTransferCoeff/liquid.Cp(), TL)
      + SpL*TV
      - fvm::Sp(SpL,TL)
      - alphaL*(fvc::ddt(rhoL,KL) + fvc::div(rhoPhiL,KL,KdivScheme))/liquid.Cp()
     //+ Sh()*mixture.rCv()
    );
    
    fvScalarMatrix TEqnV
    (
        fvm::ddt(rhoV, TV)
      + fvm::div(rhoPhiV, TV, divScheme)
      - fvm::Sp(fvc::ddt(rhoV) + fvc::div(rhoPhiV) - vapor.m_evap_sum(), TV)
      - fvm::laplacian(vapor.kByCp() + turbulence->muEff(), TV)
     ==
        heatTransferCoeff*TL/vapor.Cp()
      - fvm::Sp(heatTransferCoeff/vapor.Cp(), TV)
      + SpV*TL
      - fvm::Sp(SpV,TV)
      + alphaV*
        (
          dpdt
        - (fvc::ddt(rhoV,KV) + fvc::div(rhoPhiV,KV,KdivScheme))
        )/vapor.Cp()
     //+ Sh()*mixture.rCv()
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
