    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      + turbulence->divDevRhoReff(U)
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.interface().stf()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p)
                ) * mesh.magSf()
            )
        );
    }
    
