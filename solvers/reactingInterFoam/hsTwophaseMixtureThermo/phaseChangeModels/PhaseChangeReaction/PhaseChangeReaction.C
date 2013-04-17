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
namespace evaporationModels
{
    defineTypeNameAndDebug(PhaseChangeReaction, 0);
    addToRunTimeSelectionTable(evaporationModel, PhaseChangeReaction, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::evaporationModels::PhaseChangeReaction::PhaseChangeReaction
(
    dictionary specieDict,
    const volScalarField& p,
    const volScalarField& T,
    const phase& alphaL,
    const phase& alphaV
)
:
    evaporationModel(typeName, specieDict, p, T, alphaL, alphaV),
    reactants_(),
    products_(),
    Ta_(),
    k_()
{
}
    
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::evaporationModels::PhaseChangeReaction::calculate
(
    const volScalarField& evapMask
)
{
    evaporationModel::calculate(evapMask);

   /* m_evap_ = (coeffC_ + coeffV_) * (p_vap_*xL_ - p_*x_);

    Foam::Info<< "Min,max evaporation flux for " << vapor_specie_ << " = "
              << Foam::min(m_evap_).value() << ", " 
              << Foam::max(m_evap_).value() << " kg/m2/s" << Foam::endl;*/
}


// get the explicit and implicit source terms for Ys
Pair<tmp<volScalarField> > Foam::evaporationModels::PhaseChangeReaction::YSuSp
(
    const word& specieName
)
const
{
    tmp<volScalarField> xByY = alphaV_.xByY(vapor_specie_);
    
    //S_Yv = explicit - implicit*Yv
    // This casts the evaporation in the form C*(Ysat - Y)
        
    return Pair<tmp<volScalarField> >
    (
        area_*(coeffC_+coeffV_)*p_vap_*xL_,
        area_*(coeffC_+coeffV_)*xByY*p_
    );
}

// get the explicit and implicit source terms for pressure
Pair<tmp<volScalarField> > Foam::evaporationModels::PhaseChangeReaction::pSuSp() const
{

    return Pair<tmp<volScalarField> >
    (
        area_*(coeffC_+coeffV_)*p_vap_*xL_,
        area_*(coeffC_+coeffV_)*x_
    );
    
    /*return Pair<tmp<volScalarField> >
    (
        area_*m_evap_,
        area_*(coeffC_+coeffV_)*x_*0.0
    );*/
}


Pair<tmp<volScalarField> > Foam::evaporationModels::PhaseChangeReaction::TSuSp() const
{

    //Sh = -m_evap_*L()*area_
    //
    //Sh = explicit - implicit * T
    //
    //  implicit = -dSh/dT = (dm/dT*L() + m_evap*dL/dT)*area_
    //    dm/dT = dcoeff/dT*(pv*xL-p*x) + coeff*xL*dpv/dT
    //    dcoeff/dT = -coeff/(2*T)
    //    dpv/dT = pv * L / (R*T^2)

    tmp<volScalarField> tL = L();
    tmp<volScalarField> dmdT = -m_evap_/(2*T_) + (coeffC_+coeffV_)*xL_*dPvdT();
    tmp<volScalarField> Sp = dmdT*tL() + m_evap_*dLdT();
    tmp<volScalarField> Su = -m_evap_*tL() + Sp()*T_;
    
    return Pair<tmp<volScalarField> >
    (
        area_*Su,
        area_*Sp
    );
    
    /*return Pair<tmp<volScalarField> >
    (
        -area_*m_evap_*L(),
        -area_*m_evap_*L()/T_*0.0
    );*/
}



// ************************************************************************* //
