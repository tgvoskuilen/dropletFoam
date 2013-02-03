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

#include "evaporationModel.H"
#include "diffusionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(evaporationModel, 0);
    defineRunTimeSelectionTable(evaporationModel, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::evaporationModel::evaporationModel
(
    const word& type,
    dictionary specieDict,
    const volScalarField& p,
    const volScalarField& T,
    const phase& alphaL,
    const phase& alphaV
)
:
    alphaL_(alphaL),
    alphaV_(alphaV),
    p_(p),
    T_(T),
    evapDict_(specieDict.subDict(type + "Coeffs")),
    vapor_specie_(evapDict_.lookup("vapor")),
    m_evap_
    (
        IOobject
        (
            "m_evap_" + vapor_specie_,
            alphaL.mesh().time().timeName(),
            alphaL.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphaL.mesh(),
        dimensionedScalar("m_evap",dimDensity/dimTime,0.0)
    ),
    mask_
    (
        IOobject
        (
            "mask_" + vapor_specie_,
            alphaL.mesh().time().timeName(),
            alphaL.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphaL.mesh(),
        dimensionedScalar("mask",dimless,0.0)
    ),
    Lb_(evapDict_.lookup("Lb")),
    Tb_(evapDict_.lookup("Tb")),
    La_(0.0),
    Tc_(evapDict_.lookupOrDefault("Tc",dimensionedScalar("Tc",dimTemperature,300))),
    R_(dimensionedScalar("R", dimensionSet(1, 2, -2, -1, -1), 8314)), // J/kmol-K
    W_(alphaV_.subSpecies()[vapor_specie_]->W())
{
    if( evapDict_.found("La") )
    {
        La_ = readScalar(evapDict_.lookup("La"));
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> Foam::evaporationModel::L() const
{
    tmp<volScalarField> tL
    (
        new volScalarField
        (
            IOobject
            (
                "tL",
                alphaL_.mesh().time().timeName(),
                alphaL_.mesh()
            ),
            alphaL_.mesh(),
            Lb_ / W_
        )
    );
        
    if( La_ > 0.0 )
    {
        tL() *= Foam::pow(Foam::mag(Tc_ - T_)/(Tc_ - Tb_),La_)*pos(Tc_ - T_);
    }
    
    return tL;
}

tmp<volScalarField> Foam::evaporationModel::Sh() const
{
    return -L() * m_evap_;
}

tmp<volScalarField> Foam::evaporationModel::area() const
{
    const volScalarField& YL = alphaL_.Y(vapor_specie_+"L");
    return Foam::mag(fvc::grad(YL)) * mask_;
}


// ************************************************************************* //

