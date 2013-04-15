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
    Foam::PhaseChangeReaction

Description
    Chemical phase change

SourceFiles
    PhaseChangeReaction.C

\*---------------------------------------------------------------------------*/

#ifndef PhaseChangeReaction_H
#define PhaseChangeReaction_H

#include "evaporationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace evaporationModels
{
/*---------------------------------------------------------------------------*\
                      Class evaporationModel Declaration
\*---------------------------------------------------------------------------*/

class PhaseChangeReaction
:
    public evaporationModel
{

    //model-specific stuff

        //tables of specie names and stoich coefficients
        HashTable<scalar> reactants_;
        HashTable<scalar> products_;
        
        dimensionedScalar Ta_;
        dimensionedScalar k_;
        
        

public:

    //- Runtime type information
    TypeName("PhaseChangeReaction");

    // Constructors

        //- construct from components
        PhaseChangeReaction
        (
            dictionary specieDict,
            const volScalarField& p,
            const volScalarField& T,
            const phase& alphaL,
            const phase& alphaV
        );
        
    //- Destructor
        virtual ~PhaseChangeReaction() {}


    // Member Functions
        virtual void calculate
        (
            const volScalarField& evapMask
        );
        
        virtual Pair<tmp<volScalarField> > YSuSp(const word& specie) const;
        virtual Pair<tmp<volScalarField> > pSuSp() const;
        virtual Pair<tmp<volScalarField> > TSuSp() const;
};  


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace evaporationModel

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
