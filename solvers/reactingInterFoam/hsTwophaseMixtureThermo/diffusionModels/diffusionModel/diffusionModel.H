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
    Foam::diffusionModel

Description
    General diffusion model

SourceFiles
    diffusionModel.C

\*---------------------------------------------------------------------------*/

#ifndef diffusionModel_H
#define diffusionModel_H

#include "volFields.H"
#include "PtrDictionary.H"
#include "dictionaryEntry.H"
#include "fvCFD.H"
#include "dimensionedScalar.H"
#include "dimensionedScalarFwd.H"
#include "subSpecie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   
/*---------------------------------------------------------------------------*\
                      Class diffusionModel Declaration
\*---------------------------------------------------------------------------*/

class diffusionModel
{
    // Protected data
    protected:
        const word name_;
        
        dictionary diffusionDict_;
        
        
        diffusionModel(const diffusionModel&);
        
        void operator=(const diffusionModel&);
        
        
public:

    //- Runtime type information
    TypeName("diffusionModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            diffusionModel,
            components,
            (
                dictionary specieDict
            ),
            (specieDict)
        );


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<diffusionModel> New
        (
            dictionary specieDict
        );
        
    // Constructors
    
        diffusionModel
        (
            const word& type,
            dictionary specieDict
        );
        
        virtual ~diffusionModel() {}


    // Member Functions
    
        word name() const
        {
            return name_;
        }


        virtual tmp<volScalarField> Dij
        (
            const subSpecie& sI,
            const subSpecie& sJ
        ) const = 0;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
