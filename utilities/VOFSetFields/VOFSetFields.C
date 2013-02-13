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
    VOFSetFields

Description
    Set volume fraction fields that are not aligned with the mesh


\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "fvMesh.H"
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "droplet.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createFields.H"
    #include "createPhi.H" //needed to do boundary corrections.
    
    // Set drop species and liquid fraction
    forAllIter(PtrDictionary<droplet>, droplets, dropI)
    {
        dropI().set(alphaLiquid, U, species);
    }
    
    alphaLiquid.correctBoundaryConditions();
    
    // Set vapor fraction
    alphaVapor = 1.0 - alphaLiquid;
    alphaVapor.correctBoundaryConditions();
    
    // Set surrounding species
    forAll(species, i)
    {
        Yremaining -= species[i];
    }
    Yremaining.max(0.0);
    
    for(label i = 0; i < defaultSpecieValues.size(); ++i)
    {
        word name = defaultSpecieValues.toc()[i];
        scalar value = defaultSpecieValues[name];
        
        forAll(species, j)
        {
            if( species[j].name() == name )
            {
                species[j] += alphaVapor * Yremaining * value;
                break;
            }
        }
        
    }
    
    forAll(species, j)
    {
        species[j].max(0.0);
        species[j].correctBoundaryConditions();
    }
	
	runTime.writeNow();
	

    return 0;
}


// ************************************************************************* //
