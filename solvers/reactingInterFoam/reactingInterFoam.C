/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    ICRFoam

Description
    Solver for 2 reacting fluids, one compressible, which captures the 
    interfaces and includes surface-tension and contact-angle effects 
    for each phase.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "hsTwophaseMixtureThermo.H"
#include "turbulenceModel.H"
#include "rhoChemistryCombustionModel.H"
#include "pimpleControl.H"
#include "subCycle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "readTimeControls.H"

    scalar maxVolSource = 0.0;
    scalar MaxFo = 0.0;
    scalar volSourceLimit =
        runTime.controlDict().lookupOrDefault<scalar>("volSourceLimit", 1e-2);

    pimpleControl pimple(mesh);
    volScalarField divU(fvc::div(phi));

    #include "correctPhi.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readControls.H"
        #include "alphaCourantNo.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        //Limit deltaT based on maxVolSource
        if (adjustTimeStep)
        {
            //scalar maxVolDeltaT = maxDeltaT;

            //if( mag(maxVolSource) > SMALL )
            //    maxVolDeltaT = volSourceLimit / (mag(maxVolSource));


            if( MaxFo > 0.4 )
            {
                runTime.setDeltaT
                (
                    runTime.deltaTValue() * 0.4 / MaxFo
                );
            }
                
            //Info<< "Source*dT = " << mag(maxVolSource)*runTime.deltaTValue() 
            //    << endl;

            /*runTime.setDeltaT
            (
                min
                (
                    runTime.deltaTValue(),
                    maxVolDeltaT
                )
            );*/

            Info<< "Fo-limited deltaT = " << runTime.deltaTValue() << endl;
        }

        runTime++;
        
        // BEGIN MESH ADAPTATION SECTION
        Info<< "Time = " << runTime.timeName() << nl << endl;

        {
            // Store divU from the previous mesh for correctPhi.H
            divU = fvc::div(phi);

            //Identify regions near ANY interface or reaction zone
            refinementField = mixture.getRefinementField();
 
            // Do any mesh changes
            mesh.update();
        }

        if (mesh.changing())
        {
            gh = g & mesh.C();
            ghf = g & mesh.Cf();
        }

        if (mesh.changing() && correctPhi)
        {
            #include "correctPhi.H"
        }

        if (mesh.changing() && checkMeshCourantNo)
        {
            #include "meshCourantNo.H"
        }
        // END MESH ADAPTATION

        // This line is most likely not needed since mixture.solve
        // does this too
        // solve( fvm::ddt(rho) + fvc::div(mixture.rhoPhi()) );
        mixture.solve( rho );
        
        while (pimple.loop())
        {
            // --- Phase-Pressure-Velocity PIMPLE corrector loop
            //Info<<"Solving alpha transport equations"<<endl;
            //MaxFo = mixture.solve( rho );

            //dQ = combustion->dQ();

            #include "UEqn.H"	
            //#include "TEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }


        #include "checkMassBalance.H"
        
        //mixture.calculateRho(); //needed?
        //rho = thermo.rho();
        
        runTime.write();
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
