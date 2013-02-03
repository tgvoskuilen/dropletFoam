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
    Foam::LennardJones

Description
    LennardJones diffusion coefficient model

SourceFiles
    LennardJones.C

\*---------------------------------------------------------------------------*/

#ifndef LennardJones_H
#define LennardJones_H

#include "diffusionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusionModels
{
/*---------------------------------------------------------------------------*\
                      Class LennardJones Declaration
\*---------------------------------------------------------------------------*/

class LennardJones
:
    public diffusionModel
{

    //model-specific stuff
        scalar sigma_;
        
        scalar T_;

public:

    //- Runtime type information
    TypeName("LennardJones");

    // Constructors

        //- construct from components
        LennardJones
        (
            dictionary specieDict
        );
        
    //- Destructor
        virtual ~LennardJones() {}


    // Member Functions
        virtual tmp<volScalarField> Dij
        (
            const subSpecie& sI,
            const subSpecie& sJ
        ) const;
        
        scalar sigma() const
        {
            return sigma_;
        }
        
        scalar T() const
        {
            return T_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace diffusionModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
