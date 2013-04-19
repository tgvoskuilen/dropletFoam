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

#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseChangeModel, 0);
    defineRunTimeSelectionTable(phaseChangeModel, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeModel::phaseChangeModel
(
    const word& type,
    const word name,
    const fvMesh& mesh,
    const phase& alphaL,
    const phase& alphaV,
    const PtrList<gasThermoPhysics>& speciesData,
    dictionary phaseChangeDict
)
:
    mesh_(mesh),
    p_(mesh.lookupObject<volScalarField>("p")),
    T_(mesh.lookupObject<volScalarField>("T")),
    alphaL_(alphaL),
    alphaV_(alphaV),
    name_(name),
    phaseChangeDict_(phaseChangeDict.subDict(type + "Coeffs")),
    omega_
    (
        IOobject
        (
            "omega_"+name,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("omega_"+name, dimless/dimTime, 0.0)
    ),
    mask_
    (
        IOobject
        (
            "mask_"+name,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("mask_"+name, dimless, 0.0)
    ),
    reactants_(phaseChangeDict.lookup("liquidReactants")),
    products_(phaseChangeDict.lookup("gasProducts")),
    reacThermo_(reactants_.size()),
    prodThermo_(products_.size()),
    R_(dimensionedScalar("R", dimensionSet(1, 2, -2, -1, -1), 8314)) // J/kmol-K
{
    List<word> reacList = reactants_.toc();
    List<word> prodList = products_.toc();
    
    forAll(reacList, i)
    {
        word R = reacList[i];
        
        forAll(speciesData, s)
        {
            if( speciesData[s].name() == R )
            {
                reacThermo_.set
                (
                    R,
                    speciesData->operator()(s)
                );
                break;
            }
        }
        
        if( !reacThermo_.found(R) )
        {
            Info<< "WARNING: reacThermo for " << R << " not found in " 
                << speciesData << endl;
        }
    }
    
    forAll(prodList, i)
    {
        word P = prodList[i];
        
        forAll(speciesData, s)
        {
            if( speciesData[s].name() == P )
            {
                prodThermo_.set
                (
                    P,
                    speciesData->operator()(s)
                );
                break;
            }
        }
        
        if( !prodThermo_.found(P) )
        {
            Info<< "WARNING: prodThermo for " << P << " not found in " 
                << speciesData << endl;
        }
    }
    
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::autoPtr<Foam::phaseChangeModel> Foam::phaseChangeModel::clone() const
{
    notImplemented("phaseChangeModel::clone() const");
    return autoPtr<phaseChangeModel>(NULL);
}

bool Foam::phaseChangeModel::hasSpecie(const word& specie) const
{
    return products_.found(specie) || reactants_.found(specie);
}

Foam::tmp<Foam::volScalarField> Foam::phaseChangeModel::Sh() const
{
    Pair<tmp<volScalarField> > Ts = TSuSp();
    return Ts.first() - Ts.second()*T_;
}


// ************************************************************************* //

