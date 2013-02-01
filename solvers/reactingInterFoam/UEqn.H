
    Info<< "Building UEqn" << endl;
    
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(mixture.rhoPhi(), U)
      + turbulence->divDevRhoReff(U)
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        Info<< "Doing momentum predictor step" << endl;
        
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p)
                ) * mesh.magSf()
            )
        );
    }

