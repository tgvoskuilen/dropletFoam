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
    Foam::phaseChangeModel

Description
    General phase change model

SourceFiles
    phaseChangeModel.C

\*---------------------------------------------------------------------------*/

#ifndef phaseChangeModel_H
#define phaseChangeModel_H

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
                      Class phaseChangeModel Declaration
\*---------------------------------------------------------------------------*/

class phaseChangeModel
{
    // Private data
    protected:
        const fvMesh& mesh_;
        
        const volScalarField& p_;
        const volScalarField& T_;
        const phase& alphaL_;
        const phase& alphaV_;
        
        
        word name_;
        dictionary phaseChangeDict_;

        //- Reaction rate (kmol/m3/s) of "reactant 0"
        volScalarField omega_;
        volScalarField mask_;
        
        HashTable<scalar> reactants_;
        HashTable<scalar> products_;
        
        HashTable<const gasThermoPhysics*> reacThermo_;
        HashTable<const gasThermoPhysics*> prodThermo_;

     
        //- Universal gas constant
        dimensionedScalar R_;
        
        phaseChangeModel(const phaseChangeModel&);
        
        void operator=(const phaseChangeModel&);
        
        
public:

    //- Runtime type information
    TypeName("phaseChangeModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseChangeModel,
            components,
            (
                const word name,
                const fvMesh& mesh,
                const phase& alphaL,
                const phase& alphaV,
                const PtrList<gasThermoPhysics>& speciesData,
                dictionary phaseChangeDict
            ),
            (name, mesh, alphaL, alphaV, speciesData, phaseChangeDict)
        );


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<phaseChangeModel> New
        (
            const word name,
            const fvMesh& mesh,
            const phase& alphaL,
            const phase& alphaV,
            const PtrList<gasThermoPhysics>& speciesData,
            dictionary phaseChangeDict
        );


        //- Return clone
        autoPtr<phaseChangeModel> clone() const;
        
        //- Return a pointer to a new phase change model created on freestore
        //  from Istream
        class iNew
        {
            const fvMesh& mesh_;
            const phase& alphaL_;
            const phase& alphaV_;
            const PtrList<gasThermoPhysics>& speciesData_;

        public:

            iNew
            (
                const fvMesh& mesh,
                const phase& alphaL,
                const phase& alphaV,
                const PtrList<gasThermoPhysics>& speciesData
            )
            :
                mesh_(mesh),
                alphaL_(alphaL),
                alphaV_(alphaV),
                speciesData_(speciesData)
            {}

            autoPtr<phaseChangeModel> operator()(Istream& is) const
            {
                dictionaryEntry ent(dictionary::null, is);
                return phaseChangeModel::New
                (
                    ent.keyword(), 
                    mesh_, 
                    alphaL_, 
                    alphaV_, 
                    speciesData_
                    ent
                );
            }
        };
        
        
    // Constructors
    
        phaseChangeModel
        (
            const word& type,
            const word name,
            const fvMesh& mesh,
            const phase& alphaL,
            const phase& alphaV,
            const PtrList<gasThermoPhysics>& speciesData,
            dictionary phaseChangeDict
        );
        
        virtual ~phaseChangeModel() {}


    // Member Functions
    
        const word& name() const
        {
            return name_;
        }
        
        const word& keyword() const
        {
            return name_;
        }
        
        bool hasSpecie(const word& name) const;
        
        
    // Virtual Functions
    
        virtual void calculate
        (
            const volScalarField& phaseChangeZones,
            const volScalarField& area
        );
        
        //- Enthalpy source term
        //  Sh = TSuSp.first() - TSuSp.second()*T ?
        virtual tmp<volScalarField> Sh() const;
        
        
        
    // Abstract Functions
    

    
        //- Get the total mass generation rate for a named phase
        virtual tmp<volScalarField> mdot(const word& phaseName) const = 0;
        
        //- Get the total volume generation rate for a named phase
        virtual tmp<volScalarField> Vdot(const word& phaseName) const = 0;
               
        //- Get the implicit and explicit mass fraction source terms
        virtual Pair<tmp<volScalarField> > YSuSp(const word& name) const = 0;
        
        //- Get the implicit and explicit temperature source terms
        virtual Pair<tmp<volScalarField> > TSuSp() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
