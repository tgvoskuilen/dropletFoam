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
namespace phaseChangeModels
{
    defineTypeNameAndDebug(PhaseChangeReaction, 0);
    addToRunTimeSelectionTable(phaseChangeModel, PhaseChangeReaction, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::phaseChangeModels::PhaseChangeReaction::PhaseChangeReaction
(
    const word name,
    const fvMesh& mesh,
    const phase& alphaL,
    const phase& alphaV,
    dictionary phaseChangeDict
)
:
    phaseChangeModel(typeName, name, mesh, alphaL, alphaV, phaseChangeDict),
    Ta_(phaseChangeDict_.lookup("Ta")),
    A_(phaseChangeDict_.lookup("A")),
    beta_(readScalar(phaseChangeDict_.lookup("beta")))
{
    


    Info<< "Liquid/vapor reaction configured for " << name_ << endl;
}
    
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::phaseChangeModels::PhaseChangeReaction::calculate
(
    const volScalarField& evapMask,
    const volScalarField& area
)
{
    dimensionedScalar sA("sA",dimless/dimLength,SMALL);
    mask_ = pos(area-sA);
    
    
    omega_ = A_*Foam::exp(-Ta_/T_)*mask_;
    
    if( Foam::mag(beta_) > SMALL )
    {
        omega_ *= Foam::pow(T_,beta_);
    }
    
    dimensionedScalar c0("c0",dimMoles/dimVolume,1.0);
    List<word> R = reactants_.toc();
    forAll(R, i)
    {
        tmp<volScalarField> cL = alphaL_.x(R[i])*alphaL_.Npp()*alphaL_.rho(p_,T_)/c0;
        scalar stoic = reactants_[R[i]];
        
        if( Foam::mag(stoic-1.0) < SMALL )
        {
            omega_ *= cL;
        }
        else
        {
            omega_ *= Foam::pow(cL,stoic);
        }
    }
    
    
    Foam::Info<< "Min,max reaction flux for " << name_ << " = "
              << Foam::min(omega_).value() << ", " 
              << Foam::max(omega_).value() << " kmol/m3/s" << Foam::endl;
}


// get the explicit and implicit source terms for Ys
Pair<tmp<volScalarField> > Foam::phaseChangeModels::PhaseChangeReaction::YSuSp
(
    const word& specieName
)
const
{
    tmp<volScalarField> tSY
    (
        new volScalarField
        (
            IOobject
            (
                "SY",
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

    if( products_.found(specieName) )
    {
        tSY() += products_[specieName]*omega_*prodThermo_[specieName]->W();
    }
    else if( reactants_.found(specieName) )
    {
        tSY() -= reactants_[specieName]*omega_*reacThermo_[specieName]->W();
    }
        
    //Fully explicit. The liquid part could be made implicit, but that's not
    // the part that needs stabilizing anyway
    return Pair<tmp<volScalarField> >
    (
        tSY,
        tZero
    );
}



Pair<tmp<volScalarField> > Foam::phaseChangeModels::PhaseChangeReaction::TSuSp() const
{
    tmp<volScalarField> tSh
    (
        new volScalarField
        (
            IOobject
            (
                "Sh",
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
    
    volScalarField& Sh = tSh();
    
    forAll(tSh(), cellI)
    {
        if( mag(omega_[cellI]) > SMALL )
        {
            forAllConstIter(HashTable<const gasThermoPhysics*>, reacThermo_, tRi)
            {
                word specie = tRi()->name();
                
                const scalar hi = tRi()->hc(); //Hc = H(Tstd) in J/kmol
                
                Sh[cellI] += hi*omega_[cellI]*reactants_[specie];
            }
            
            forAllConstIter(HashTable<const gasThermoPhysics*>, prodThermo_, tPi)
            {
                word specie = tPi()->name();
                
                const scalar hi = tPi()->hc(); //Hc = H(Tstd) in J/kmol
                
                Sh[cellI] -= hi*omega_[cellI]*products_[specie];
            }
        }
    }

    return Pair<tmp<volScalarField> >
    (
        tSh,
        tZero
    );
}

//net mass generation in phase 'phaseName'
tmp<volScalarField> Foam::phaseChangeModels::PhaseChangeReaction::mdot
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

tmp<volScalarField> Foam::phaseChangeModels::PhaseChangeReaction::Vdot
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
