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
    reactingInterFoam

Description
    Solver for 2 reacting fluids, one compressible, which captures the 
    interfaces and includes surface-tension.

    Turbulence modeling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "hsTwophaseMixtureThermo.H"
#include "turbulenceModel.H"
#include "rhoChemistryCombustionModel.H"
#include "pimpleControl.H"
#include "subCycle.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "readTimeControls.H"

    scalar totalMass = 0.0;
    scalar massError = 0.0;

    word clipDir = runTime.controlDict().lookupOrDefault<word>
    (
        "clipDirection",
        "none"
    );

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
            mixture.updateMeshArDelta();

            if (correctPhi)
            {
                #include "correctPhi.H"
            }

            if (checkMeshCourantNo)
            {
                #include "meshCourantNo.H"
            }
        }
        // END MESH ADAPTATION
        
        while (pimple.loop())
        {
            // --- Phase-Pressure-Velocity PIMPLE corrector loop
            Info<<"Solving alpha transport equations"<<endl;
            mixture.solve( rho, pimple.corr() );

            dQ = combustion->dQ() + mixture.dQ_phaseChange();

            #include "UEqn.H"	
            #include "TEqn.H"
            
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
        #include "writeSummaryParameters.H"
        
        runTime.write();
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
