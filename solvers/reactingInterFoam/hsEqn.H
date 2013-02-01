{   
    // TODO: Change this to TEqn?
    //   ddt term - get rhoCp from mixture
    //   convection - explicit convection from adding phases (or implicit with Cp)
    //   laplacian - implicit, use mixture.k()
    //   source terms - should require no change, W/m3 is correct unit

    // then it can be clipped at 200-6000 K before calculating hs


    // TODO: Solve hs for each phase?
    
    dimensionedScalar keff("keff",dimPower/dimLength/dimTemperature, 0.1);
    
    //volScalarField rhoNew(thermo.rho());
    //solve(fvm::ddt(rhoNew) + fvc::div(mixture.rhoPhi()) - divU*rho);
    
    fvScalarMatrix hsEqn
    (
        fvm::ddt(rho, hs)
      + fvc::div(mixture.rhoPhi(), hs) //TODO: Have each phase store its rhoPhi and calculate rhoPhiHs
      //- divU*rhoNew*hs
      //+ fvc::div(mixture.rhoPhiHs())
      //- fvm::laplacian(mixture.alpha(), hs) //Diffusion of hs is non-physical in a multiphase system...
      //- fvm::laplacian(turbulence->alphaEff(), hs)
      //- fvc::laplacian(mixture.k_Cp(), mixture.Hc())
      - fvc::laplacian(keff, T)
     ==
        dpdt
      - (fvc::ddt(rho, K) + fvc::div(mixture.rhoPhi(), K))
      + Sh()
    );

    hsEqn.relax();
    hsEqn.solve();

    // Update temperature field using new hs
    thermo.correct();

    Info<< "T min/max   = "
        << gMin(T) << ", "
        << gMax(T) << endl;

}

