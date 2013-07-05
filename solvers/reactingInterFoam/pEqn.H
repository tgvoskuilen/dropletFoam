{
    const volScalarField& alphaV = mixture.alphaVapor();
    dimensionedScalar pMin("pMin",p.dimensions(),10000.0);
    
    rho = thermo.rho();
    
    // First part of density update, only changes compressible regions
    //  since psi is psiV*alphaV and is zero in the liquid regions
    thermo.rho() -= thermo.psi() * p_rgh;
    
    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rAUf(fvc::interpolate(rAU));

    // Dp/Dt term
    //   compressibleInterFoam uses DrhoDt here, fully explicit, rather than
    //   separate out the compressibility factor and make time-implicit.
    //   Transonic cases make this fully implicit.
    /*fvScalarMatrix p_DpDt = 
    (
        fvm::ddt(p)
      + fvc::div(phi, p)
      - fvc::Sp(fvc::div(phi), p)
    );*/
    
   
    // Estimate U with no force terms
    U = rAU*UEqn.H();
    
    if(clipDir=="x")
    {
        forAll(U,cellI)
        {
            U[cellI].x() = 0.0;
        }
    }
    else if(clipDir=="y")
    {
        forAll(U,cellI)
        {
            U[cellI].y() = 0.0;
        }
    }
    else if(clipDir=="z")
    {
        forAll(U,cellI)
        {
            U[cellI].z() = 0.0;
        }
    }
    
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(U) & mesh.Sf())
      + fvc::ddtPhiCorr(rAU, rho, U, phi)
    );

    // Add surface tension and gravity forces to flux estimate
    //phi = phiHbyA +
    surfaceScalarField phiForces
    (
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );
    
    phi = phiHbyA + phiForces;
    
    // DpDt term, uses phiEst based on buoyantPimpleFoam pEqn.H
    fvScalarMatrix p_rghDDtEqn = 
    (
        fvm::ddt(p_rgh)
      + fvc::div(phi, p_rgh)
      - fvc::Sp(fvc::div(phi), p_rgh)
    );
    
    tmp<volScalarField> p_rgh_min = pMin - p0 - rho*gh;
     
    // Solve for pressure to enforce continuity (based on phi)
    while (pimple.correctNonOrthogonal())
    {   
        fvScalarMatrix p_rghEqn
        (
            fvc::div(phi)
          - fvm::laplacian(rAUf, p_rgh)
          - mixture.divPhaseChange()
        );
        
        solve
        (
            max(alphaV,0.0)/p*p_rghDDtEqn
          - max(alphaV,0.0)/T*fvc::DDt(phi,T)
          + p_rghEqn,
          
            mesh.solver(p_rgh.select(pimple.finalInnerIter()))
        );
            
        p = p_rgh + p0 + rho*gh;
        p.max(pMin);
            
        if (pimple.finalNonOrthogonalIter())
        {
            // Second part of thermodynamic density update
            thermo.rho() += thermo.psi()*p_rgh;

            phi += p_rghEqn.flux();
            //p_rgh.relax(); //from Euler solver and Buoyant solver...
            
            U += rAU*fvc::reconstruct((phiForces + p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions(); 
    
            //should be fvc::DDt(phi,psi)/psi but psi is 0 some places
            mixture.divComp() = (p_rghDDtEqn & p_rgh)/p
                              - fvc::DDt(phi,T)/T;
        }
    }
    
   
    
    p = p_rgh + p0 + rho*gh;
    p.max(pMin);
    
    K = 0.5*magSqr(U);
    dpdt = fvc::ddt(p_rgh);
    
    // Continuity error check checks whether thermo.rho satisfies continuity
    // with rhoPhi field
    // However, since we aren't updated rhoPhi, this would not make sense
    // so we don't check, and we instead just call thermo.correct
    //thermo.correct();
    
    Info<< "min,max p = " << gMin(p)/1e5 << ", " 
        << gMax(p)/1e5 << " bar" << endl;
        
    Info<< "min,max U = " << min(mag(U)).value() 
        << ", " << max(mag(U)).value() << endl;
}
