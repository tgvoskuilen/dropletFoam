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
    Foam::evaporationModel

Description
    General evaporation model

SourceFiles
    evaporationModel.C

\*---------------------------------------------------------------------------*/

#ifndef evaporationModel_H
#define evaporationModel_H

#include "volFields.H"
#include "PtrDictionary.H"
#include "dictionaryEntry.H"
#include "fvCFD.H"
#include "dimensionedScalar.H"
#include "dimensionedScalarFwd.H"
#include "phase.H"
#include "subSpecie.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   
/*---------------------------------------------------------------------------*\
                      Class evaporationModel Declaration
\*---------------------------------------------------------------------------*/

class evaporationModel
{
    // Private data
    protected:
        const phase& alphaL_;
        const phase& alphaV_;
        const volScalarField& p_;
        const volScalarField& T_;
        dictionary evapDict_;

        word vapor_specie_;

        volScalarField m_evap_;
        //volScalarField mask_;
        
        //- Heat of vaporization (J/kmol)
        //  L = Lb * ((Tc-T)/(Tc-Tb))^La
        dimensionedScalar Lb_;
        dimensionedScalar Tb_;
        scalar La_;

        //- Specie critical temperature
        dimensionedScalar Tc_;
                
        //- Universal gas constant
        dimensionedScalar R_;
        
        //- Specie molar mass
        dimensionedScalar W_;
        
        
        evaporationModel(const evaporationModel&);
        
        void operator=(const evaporationModel&);
        
        
public:

    //- Runtime type information
    TypeName("evaporationModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            evaporationModel,
            components,
            (
                dictionary specieDict,
                const volScalarField& p,
                const volScalarField& T,
                const phase& alphaL,
                const phase& alphaV
            ),
            (specieDict, p, T, alphaL, alphaV)
        );


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<evaporationModel> New
        (
            dictionary specieDict,
            const volScalarField& p,
            const volScalarField& T,
            const phase& alphaL,
            const phase& alphaV
        );
        
    // Constructors
    
        evaporationModel
        (
            const word& type,
            dictionary specieDict,
            const volScalarField& p,
            const volScalarField& T,
            const phase& alphaL,
            const phase& alphaV
        );
        
        virtual ~evaporationModel() {}


    // Member Functions
    
        word name() const
        {
            return vapor_specie_+" Evaporation";
        }

        const word& vaporName() const
        {
            return vapor_specie_;
        }
        
        const volScalarField& m_evap() const
        {
            return m_evap_;
        }
        
        virtual void calculate
        (
            const volScalarField& area
        ) = 0;
        
        //- Calculate enthalpy source due to evaporation
        tmp<volScalarField> Sh() const;
        
        //- Calculate latent heat of vaporization (function of T)
        tmp<volScalarField> L() const;
        
        //- Calculate interface area/volume
        //tmp<volScalarField> area() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
