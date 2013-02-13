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
    multiphaseInterFoam

Description
    Solver for n incompressible fluids which captures the interfaces and
    includes surface-tension and contact-angle effects for each phase.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "fvMesh.H"
#include "fvCFD.H"
#include "dynamicFvMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	
    Info<< "Reading field alpha1\n" << endl;

    /*
        The mesh refinement uses 'refinementField' as the refinement criteria.
        In my solvers, 'refinementField' is calculated from the multiphase
        mixture.

        In this pre-adaptation routine, 'refinementField' is calculated from
        the interfaces of alphaair, as specified in setFields.

        Use this routine iteratively, alternating with setFields to refine and
        capture the initial conditions
            initDynamicMesh
            copy originals from 0.org to new time folder
            setFields -latestTime


    */

    volScalarField alphaVapor
    (
        IOobject
        (
            "alphaVapor",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    //relax alphaVapor a little bit
    /*dimensionedScalar dx = pow(min(mesh.V()), 1.0/3.0);
    dimensionedScalar dA = dx*dx/4.0; //dTau("dTau",dimArea,dx*dx/4.0);
    
    for(label i = 0; i < 1; i++)
    {
        alphaVapor += dA * fvc::laplacian(alphaVapor);
    }*/


    volScalarField refinementField
    (
        IOobject
        (
            "refinementField",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        1e6*mag(fvc::grad(alphaVapor))
    );
    //refinementField.internalField() *= pow(mesh.V(),1.0/3.0);

    runTime.setDeltaT(1e-6);
    runTime++;
    
	Info<< "Time = " << runTime.timeName() << nl << endl;
	

	scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();
	{
		mesh.update();

        refinementField = 1e6*mag(fvc::grad(alphaVapor));
        mesh.update();
	}
	
	

	if (mesh.changing())
	{
		Info<< "Execution time for mesh.update() = "
			<< runTime.elapsedCpuTime() - timeBeforeMeshUpdate
			<< " s" << endl;
	}

	
	runTime.writeNow();


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
