{
    //word alphaScheme("div(phi,alpha)");
    //word alpharScheme("div(phirb,alpha)");

    //surfaceScalarField phic(mag(phi/mesh.magSf()));
    //phic = min(interface.cAlpha()*phic, max(phic));
    //surfaceScalarField phir(phic*interface.nHatf());

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        #include "calculateInterfacePlanes.H"

        #include "setPhiAlpha.H"
        surfaceScalarField phiAlpha = phiAlphaFrac * phi;

        Info<< "phiAlpha min,max = " << Foam::min(phiAlpha) << ", " << Foam::max(phiAlpha) << endl;
        phiAlphaDiv = fvc::div(phiAlpha) * runTime.deltaT();

        MULES::explicitSolve
        (
            geometricOneField(),
            alpha1, 
            phi, 
            phiAlpha, 
            zeroField(), //Sp
            zeroField(), //Su
            1,           //alphaMax
            0            //alphaMin
        );
        alpha1.max(0.0);
        Info<< "phiAlpha after MULES min,max = " << Foam::min(phiAlpha) << ", " << Foam::max(phiAlpha) << endl;

        rhoPhi = phiAlpha*(rho1 - rho2) + phi*rho2;
    }

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;
}
