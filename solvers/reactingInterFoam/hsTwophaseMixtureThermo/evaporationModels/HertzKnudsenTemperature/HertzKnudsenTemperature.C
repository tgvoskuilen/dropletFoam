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

#include "HertzKnudsenTemperature.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace evaporationModels
{
    defineTypeNameAndDebug(HertzKnudsenTemperature, 0);
    addToRunTimeSelectionTable(evaporationModel, HertzKnudsenTemperature, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::evaporationModels::HertzKnudsenTemperature::HertzKnudsenTemperature
(
    dictionary specieDict,
    const volScalarField& p,
    const volScalarField& T,
    const phase& alphaL,
    const phase& alphaV
)
:
    evaporationModel(typeName, specieDict, p, T, alphaL, alphaV),
    betaM_(readScalar(evapDict_.lookup("betaM")))
{
}
    
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::evaporationModels::HertzKnudsenTemperature::calculate
(
    const volScalarField& evapMask
)
{
    setMask( evapMask );
        
    scalar pi = constant::mathematical::pi;
    tmp<volScalarField> rhoV = p_*W_/(R_*T_);
    
    m_evap_ = area() * 2.0*betaM_/(2.0-betaM_)*Foam::sqrt(W_/(2*pi*R_))
              * L()*rhoV*(T_ - Tb_)/Foam::pow(Tb_,1.5);
    m_evap_.max(0.0);
    
    Foam::Info << "Min,max evaporation source for " << vapor_specie_ << " = "
               << Foam::gMin(m_evap_) << ", " 
               << Foam::gMax(m_evap_) << " kg/m3/s" << Foam::endl;
}

// ************************************************************************* //
