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

Class
    Foam::UnitySchmidtNo

Description
    UnitySchmidtNo diffusion coefficient model

SourceFiles
    UnitySchmidtNo.C

\*---------------------------------------------------------------------------*/

#ifndef UnitySchmidtNo_H
#define UnitySchmidtNo_H

#include "diffusionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusionModels
{
/*---------------------------------------------------------------------------*\
                      Class UnitySchmidtNo Declaration
\*---------------------------------------------------------------------------*/

class UnitySchmidtNo
:
    public diffusionModel
{

    //model-specific stuff
        

public:

    //- Runtime type information
    TypeName("UnitySchmidtNo");

    // Constructors

        //- construct from components
        UnitySchmidtNo
        (
            dictionary specieDict
        );
        
    //- Destructor
        virtual ~UnitySchmidtNo() {}


    // Member Functions
        virtual tmp<volScalarField> Dij
        (
            const subSpecie& sI,
            const subSpecie& sJ,
            const compressible::turbulenceModel& turb
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace diffusionModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
