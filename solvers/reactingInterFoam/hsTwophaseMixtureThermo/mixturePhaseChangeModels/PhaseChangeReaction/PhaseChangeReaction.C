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

\*---------------------------------------------------------------------------*/

#include "PhaseChangeReaction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mixturePhaseChangeModels
{
    defineTypeNameAndDebug(PhaseChangeReaction, 0);
    addToRunTimeSelectionTable(mixturePhaseChangeModel, PhaseChangeReaction, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::mixturePhaseChangeModels::PhaseChangeReaction::PhaseChangeReaction
(
    const word name,
    const fvMesh& mesh,
    const phase& alphaL,
    const phase& alphaV,
    const PtrList<gasThermoPhysics>& speciesData,
    dictionary phaseChangeDict
)
:
    mixturePhaseChangeModel
    (
        typeName, 
        name, 
        mesh, 
        alphaL, 
        alphaV, 
        speciesData, 
        phaseChangeDict
    ),
    key_specie_(phaseChangeDict_.lookup("keySpecie")),
    keyID_(alphaL_.subSpecies()[key_specie_]->idx()),
    keyW_(alphaL_.subSpecies()[key_specie_]->W())
{
    Info<< "Liquid/vapor reaction configured for " << name_ << endl;
}
    
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //       
void Foam::mixturePhaseChangeModels::PhaseChangeReaction::calculate
(
    const volScalarField& evapMask,
    const volScalarField& area
)
{
    //get the reaction rate based on the total mass rate of change of the
    // key specie. This assumes the key specie is only involved in this
    // reaction!
    
    const rhoChemistryModel& chemistry = combustionPtr_->pChemistry();
    tmp<volScalarField> omegaKey = chemistry.RR( keyID_ ) / keyW_;

    omega_ = -omegaKey / reactants_[key_specie_];

    Foam::Info<< "Min,max reaction flux for " << name_ << " = "
              << Foam::min(omega_).value() << ", " 
              << Foam::max(omega_).value() << " kmol/m3/s" << Foam::endl;
}


// get the explicit and implicit source terms for Ys
Pair<tmp<volScalarField> > Foam::mixturePhaseChangeModels::PhaseChangeReaction::YSuSp
(
    const word& specieName
)
const
{
    tmp<volScalarField> tZero
    (
        new volScalarField
        (
            IOobject
            (
                "zero",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimDensity/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    // Y source is already counted in reactions, so don't count it here
    return Pair<tmp<volScalarField> >
    (
        tZero,
        tZero
    );
}



Pair<tmp<volScalarField> > Foam::mixturePhaseChangeModels::PhaseChangeReaction::TSuSp() const
{    
    tmp<volScalarField> tZero
    (
        new volScalarField
        (
            IOobject
            (
                "zero",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    // Heat source is already counted in reactions, so don't count it here
    return Pair<tmp<volScalarField> >
    (
        tZero,
        tZero
    );
}

//net mass generation in phase 'phaseName'
tmp<volScalarField> Foam::mixturePhaseChangeModels::PhaseChangeReaction::mdot
(
    const word& phaseName
) const
{
    tmp<volScalarField> tmdot
    (
        new volScalarField
        (
            IOobject
            (
                "tmdot",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("tmdot",dimDensity/dimTime,0.0)
        )
    );
    
    if( phaseName == "Vapor" )
    {
        forAllConstIter(HashTable<const gasThermoPhysics*>, prodThermo_, tPi)
        {
            word specie = tPi()->name();
            tmdot() += products_[specie]*omega_*tPi()->W();
        }
    }
    else if( phaseName == "Liquid" )
    {
        forAllConstIter(HashTable<const gasThermoPhysics*>, reacThermo_, tRi)
        {
            word specie = tRi()->name();
            tmdot() -= reactants_[specie]*omega_*tRi()->W();
        }
    }
    
    return tmdot;
}

tmp<volScalarField> Foam::mixturePhaseChangeModels::PhaseChangeReaction::Vdot
(
    const word& phaseName
) const
{
    tmp<volScalarField> tVdot
    (
        new volScalarField
        (
            IOobject
            (
                "tVdot",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("tVdot",dimless/dimTime,0.0)
        )
    );
    
    if( phaseName == "Vapor" )
    {
        forAllConstIter(HashTable<const gasThermoPhysics*>, prodThermo_, tPi)
        {
            word specie = tPi()->name();
            tVdot() += products_[specie]*omega_*R_*T_/p_;
        }
    }
    else if( phaseName == "Liquid" )
    {
        forAllConstIter(HashTable<const gasThermoPhysics*>, reacThermo_, tRi)
        {
            word specie = tRi()->name();
            dimensionedScalar rhoL = alphaL_.subSpecies()[specie]->rho0();
            tVdot() -= reactants_[specie]*omega_*tRi()->W()/rhoL;
        }
    }
    
    return tVdot;
}
// ************************************************************************* //
