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
    const PtrList<gasHThermoPhysics>& speciesData,
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
    )
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
    //do nothing
}


// get the explicit and implicit source terms for Ys
Pair<tmp<volScalarField> > 
Foam::mixturePhaseChangeModels::PhaseChangeReaction::YSuSp
(
    const word& specieName
) const
{
    tmp<volScalarField> tZero1
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
    
    tmp<volScalarField> tZero2
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
        tZero1,
        tZero2
    );
}



Pair<tmp<volScalarField> > 
Foam::mixturePhaseChangeModels::PhaseChangeReaction::TSuSp() const
{    
    tmp<volScalarField> tZero1
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
            dimensionedScalar("zero", dimPower/dimVolume, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    tmp<volScalarField> tZero2
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
            dimensionedScalar("zero", dimPower/dimVolume/dimTemperature, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    // Heat source is already counted in reactions, so don't count it here
    return Pair<tmp<volScalarField> >
    (
        tZero1,
        tZero2
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
    
    if(combustionPtr_)
    {
        const rhoChemistryModel& chemistry = combustionPtr_->chemistry();
        
        if( phaseName == "Vapor" )
        {
            forAllConstIter(PtrDictionary<subSpecie>, alphaV_.subSpecies(), specieI)
            {
                tmdot().dimensionedInternalField() += chemistry.RR(specieI().idx());
            }
        }
        else if( phaseName == "Liquid" )
        {
            forAllConstIter(PtrDictionary<subSpecie>, alphaL_.subSpecies(), specieI)
            {
                tmdot().dimensionedInternalField() += chemistry.RR(specieI().idx());
            }
        }
    }
    return tmdot;
}

tmp<volScalarField> Foam::mixturePhaseChangeModels::PhaseChangeReaction::Vdot
(
    const word& phaseName
) const
{
    // While the mass divergence is 0 within a phase, the velocity divergence
    // is not guaranteed to be so, since reactions change the number of moles
    // present which changes the partial volumes of the species involved. This
    // term can be nonzero even when there are no phase-change reactions.

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
    
    if(combustionPtr_)
    {
        const rhoChemistryModel& chemistry = combustionPtr_->chemistry();
        
        if( phaseName == "Vapor" )
        {
            forAllConstIter(PtrDictionary<subSpecie>, alphaV_.subSpecies(), specieI)
            {
                tmp<volScalarField> rhoV = p_*specieI().W()/(R_*T_);
                tVdot().dimensionedInternalField() += chemistry.RR(specieI().idx())/rhoV().dimensionedInternalField();
            }
        }
        else if( phaseName == "Liquid" )
        {
            forAllConstIter(PtrDictionary<subSpecie>, alphaL_.subSpecies(), specieI)
            {
                tVdot().dimensionedInternalField() += chemistry.RR(specieI().idx())/specieI().rho0();
            }
        }
    }
    return tVdot;
}
// ************************************************************************* //
