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
    Foam::mixturePhaseChangeModel

Description
    General phase change model

SourceFiles
    mixturePhaseChangeModel.C

\*---------------------------------------------------------------------------*/

#ifndef mixturePhaseChangeModel_H
#define mixturePhaseChangeModel_H

#include "volFields.H"
#include "PtrDictionary.H"
#include "dictionaryEntry.H"
#include "fvCFD.H"
#include "dimensionedScalar.H"
#include "dimensionedScalarFwd.H"
#include "phase.H"
#include "subSpecie.H"
#include "rhoChemistryCombustion.H"
#include "rhoCombustionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   
/*---------------------------------------------------------------------------*\
                      Class mixturePhaseChangeModel Declaration
\*---------------------------------------------------------------------------*/

class mixturePhaseChangeModel
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
        
        HashTable<const gasHThermoPhysics*> reacThermo_;
        HashTable<const gasHThermoPhysics*> prodThermo_;

        //- Pointer to the combustion model (given after construction)
        combustionModels::rhoChemistryCombustion* combustionPtr_;
        
        //- Universal gas constant
        dimensionedScalar R_;
        
        mixturePhaseChangeModel(const mixturePhaseChangeModel&);
        
        void operator=(const mixturePhaseChangeModel&);
        
        
public:

    //- Runtime type information
    TypeName("mixturePhaseChangeModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            mixturePhaseChangeModel,
            components,
            (
                const word name,
                const fvMesh& mesh,
                const phase& alphaL,
                const phase& alphaV,
                const PtrList<gasHThermoPhysics>& speciesData,
                dictionary phaseChangeDict
            ),
            (name, mesh, alphaL, alphaV, speciesData, phaseChangeDict)
        );


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<mixturePhaseChangeModel> New
        (
            const word name,
            const fvMesh& mesh,
            const phase& alphaL,
            const phase& alphaV,
            const PtrList<gasHThermoPhysics>& speciesData,
            dictionary phaseChangeDict
        );


        //- Return clone
        autoPtr<mixturePhaseChangeModel> clone() const;
        
        //- Return a pointer to a new phase change model created on freestore
        // from Istream
        class iNew
        {
            const fvMesh& mesh_;
            const phase& alphaL_;
            const phase& alphaV_;
            const PtrList<gasHThermoPhysics>& speciesData_;

        public:

            iNew
            (
                const fvMesh& mesh,
                const phase& alphaL,
                const phase& alphaV,
                const PtrList<gasHThermoPhysics>& speciesData
            )
            :
                mesh_(mesh),
                alphaL_(alphaL),
                alphaV_(alphaV),
                speciesData_(speciesData)
            {}

            autoPtr<mixturePhaseChangeModel> operator()(Istream& is) const
            {
                dictionaryEntry ent(dictionary::null, is);
                return mixturePhaseChangeModel::New
                (
                    ent.keyword(),
                    mesh_,
                    alphaL_,
                    alphaV_,
                    speciesData_,
                    ent
                );
            }
        };
        
    // Constructors
    
        mixturePhaseChangeModel
        (
            const word& type,
            const word name,
            const fvMesh& mesh,
            const phase& alphaL,
            const phase& alphaV,
            const PtrList<gasHThermoPhysics>& speciesData,
            dictionary phaseChangeDict
        );
        
        virtual ~mixturePhaseChangeModel() {}


    // Member Functions
        // defined here, cannot be overwritten
        
        const word& name() const
        {
            return name_;
        }
        
        const word& keyword() const
        {
            return name_;
        }
        
        bool hasSpecie(const word& name) const;
        
        void setPtr(combustionModels::rhoCombustionModel* combustion)
        {
            combustionPtr_ = dynamic_cast<combustionModels::rhoChemistryCombustion*>(combustion);
        }
        
    // Virtual Functions
        // defined here, can be overwritten

        
        //- Enthalpy source term
        //  Sh = TSuSp.first() - TSuSp.second()*T ?
        virtual tmp<volScalarField> Sh() const;
        
        
        
    // Abstract Functions
        // not defined here, must be defined by derived class
        
        virtual void calculate
        (
            const volScalarField& phaseChangeZones,
            const volScalarField& area
        ) =0;
    
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
