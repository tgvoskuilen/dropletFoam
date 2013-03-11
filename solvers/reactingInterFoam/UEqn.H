
    Info<< "Building UEqns" << endl;
    
    word divSchemeU("div(rho*phi,U)");

    fvVectorMatrix UEqnL
    (
        fvm::ddt(rhoL, UL)
      + fvm::div(rhoPhiL, UL, divSchemeU)
      - fvm::Sp(fvc::ddt(rhoL) + fvc::div(rhoPhiL) - liquid.m_evap_sum(), UL)
      - fvm::laplacian(liquid.mu(p), UL)
      //+ fvm::Sp(SpL, UL)
      //+ turbulence->divDevRhoReff(U)
    );

    fvVectorMatrix UEqnV
    (
        fvm::ddt(rhoV, UV)
      + fvm::div(rhoPhiV, UV, divSchemeU)
      - fvm::Sp(fvc::ddt(rhoV) + fvc::div(rhoPhiV) - vapor.m_evap_sum(), UV)
      - fvm::laplacian(vapor.mu(p), UV)
      //+ fvm::Sp(SpV, UV)
      //+ turbulence->divDevRhoReff(UV)
    );
    
    UEqnL.relax();
    UEqnV.relax();

