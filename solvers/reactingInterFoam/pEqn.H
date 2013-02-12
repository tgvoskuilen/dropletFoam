{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("Dp", fvc::interpolate(rAU));

    volVectorField Uest("Uest", U);
    Uest = rAU*UEqn.H();
    
    forAll(Uest, cellI)
    {
        Uest[cellI].x() = 0.0;
    }

    surfaceScalarField phiEst
    (
        "phiEst",
        (fvc::interpolate(Uest) & mesh.Sf())
      + fvc::ddtPhiCorr(rAU, rho, U, phi)
    );

    //adjustPhi(phiHbyA, U, p); //?
    phi = phiEst;

    surfaceScalarField phig
    (
        (
            mixture.interface().stf()
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiEst += phig;

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAUf, p) == fvc::div(phiEst)
        );

        //pEqn.setReference(pRefCell, getRefCellValue(p, pRefCell));

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiEst - pEqn.flux();
            
            volScalarField Uestperp = Uest & vector(1,0,0);
            volScalarField Ucorrperp = (rAU*fvc::reconstruct((phig - pEqn.flux())/rAUf)) & vector(1,0,0);
            
            if( Foam::max(Foam::mag(Uestperp)).value() > 1 || Foam::max(Foam::mag(Ucorrperp)).value() > 1 )
            {
                Info<<"Uestperp = " << Foam::max(Foam::mag(Uestperp)) << endl;
                Info<<"Ucorrperp = " << Foam::max(Foam::mag(Ucorrperp)) << endl;
            }
            
            U = Uest + rAU*fvc::reconstruct((phig - pEqn.flux())/rAUf);
            U.correctBoundaryConditions();
        }
    }
    
    Info<<"max U = " << Foam::max(Foam::mag(U)).value() << endl;
    Info<<"max p = " << Foam::max(p).value() << endl;
    Info<<"max phi = " << Foam::max(Foam::mag(phi)).value() << endl;
    
    #include "continuityErrs.H"

    /*p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }*/
}
