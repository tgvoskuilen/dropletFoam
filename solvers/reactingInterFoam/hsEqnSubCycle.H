{
    label nHsSubCycles(readLabel(pimple.dict().lookup("nHsSubCycles")));

    tmp<volScalarField> Sh = mixture.Sh(); //combustion->Sh() + mixture.Sh_evap();
    volScalarField divU(fvc::div(phi));
    
    //thermo.correct();

    Info<< "Min,max net hs source = "
        << Foam::gMin(Sh()) << ", " << Foam::gMax(Sh()) << endl;

    if (nHsSubCycles > 1)
    {
        dimensionedScalar totalDeltaT = runTime.deltaT();

        for
        (
            subCycle<volScalarField> hsSubCycle(hs, nHsSubCycles);
            !(++hsSubCycle).end();
        )
        {
            #include "hsEqn.H"
        }
    }
    else
    {
        #include "hsEqn.H"
    }

    hsE = hs;
}


