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

#include "HertzKnudsenPressure.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace evaporationModels
{
    defineTypeNameAndDebug(HertzKnudsenPressure, 0);
    addToRunTimeSelectionTable(evaporationModel, HertzKnudsenPressure, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::evaporationModels::HertzKnudsenPressure::HertzKnudsenPressure
(
    dictionary specieDict,
    const volScalarField& p,
    const volScalarField& T,
    const phase& alphaL,
    const phase& alphaV
)
:
    evaporationModel(typeName, specieDict, p, T, alphaL, alphaV),
    x_
    (
        IOobject
        (
            "x_" + vapor_specie_,
            alphaL.mesh().time().timeName(),
            alphaL.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphaL.mesh(),
        dimensionedScalar("x",dimless,0.0)
    ),
    p_vap_
    (
        IOobject
        (
            "p_vap_" + vapor_specie_,
            alphaL.mesh().time().timeName(),
            alphaL.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphaL.mesh(),
        dimensionedScalar("pv",dimPressure,0.0)
    ),
    Pc_(evapDict_.lookup("Pc")),
    PvCoeffs_(evapDict_.lookup("PvCoeffs")),
    betaM_(readScalar(evapDict_.lookup("betaM")))
{
}
    
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::evaporationModels::HertzKnudsenPressure::calculate
(
    const volScalarField& evapMask
)
{
    //Calculate the vapor pressure
    dimensionedScalar p0("p0",dimPressure,101325);
    
    p_vap_ = p0*exp
     (
        (-Lb_/R_*(1.0/T_ - 1.0/Tb_))*neg(T_ - Tb_)

      + (
           PvCoeffs_[0]
         - dimensionedScalar("B",dimTemperature,PvCoeffs_[1])/T_
         - dimensionedScalar("C",dimTemperature*dimTemperature,PvCoeffs_[2])
         / (T_*T_)
        )*pos(T_ - Tb_)
     );
    
    p_vap_.min(Pc_);
    
    //Update the mask
    mask_ = (pos(0.8 - alphaL_) * evapMask * neg(T_ - Tc_*0.9)
             + pos(0.8 - alphaL_) * pos(T_ - Tc_*0.9)) 
             * pos(alphaV_.Yp() - 1e-3);
                          
                          
    //Calculate the vapor mol fraction
    x_ = alphaV_.Y(vapor_specie_) / (W_ * alphaV_.Np()) * mask_;
    
    scalar pi = constant::mathematical::pi;
    
    m_evap_ = area() * 2.0*betaM_/(2.0-betaM_)*Foam::sqrt(W_/(2*pi*R_*T_))
                * (p_vap_ - p_*x_);
                
    m_evap_.max(0.0);
    
    Foam::Info << "Min,max evaporation source for " << vapor_specie_ << " = "
               << Foam::gMin(m_evap_) << ", " 
               << Foam::gMax(m_evap_) << " kg/m3/s" << Foam::endl;
}

// ************************************************************************* //
