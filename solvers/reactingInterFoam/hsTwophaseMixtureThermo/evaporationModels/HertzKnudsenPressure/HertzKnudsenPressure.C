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
    xL_
    (
        IOobject
        (
            "xL_" + vapor_specie_,
            alphaL.mesh().time().timeName(),
            alphaL.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphaL.mesh(),
        dimensionedScalar("xL",dimless,0.0)
    ),
    coeffC_
    (
        IOobject
        (
            "coeffC_" + vapor_specie_,
            alphaL.mesh().time().timeName(),
            alphaL.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphaL.mesh(),
        dimensionedScalar("coeffC",dimTime/dimLength,0.0)
    ),
    coeffV_
    (
        IOobject
        (
            "coeffV_" + vapor_specie_,
            alphaL.mesh().time().timeName(),
            alphaL.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphaL.mesh(),
        dimensionedScalar("coeffV",dimTime/dimLength,0.0)
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
    evaporationModel::calculate(evapMask);
    
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
                              
    //Calculate the mole fractions

    xL_ = alphaL_.x(vapor_specie_+"L");
    
    tmp<volScalarField> xSat = p_vap_/p_*xL_;
    xSat().min(1.0);
    x_ = alphaV_.x(vapor_specie_)*pos(alphaV_.Yp() - 1e-4) + xSat*neg(alphaV_.Yp() - 1e-4);
    
    //x_ = alphaV_.x(vapor_specie_);
    
    scalar pi = constant::mathematical::pi;
    dimensionedScalar sA("sA",dimless/dimLength,SMALL);
    
    tmp<volScalarField> coeff = Foam::sqrt(W_/(2*pi*R_*T_))*pos(area_-sA)*mask_;
    
    dimensionedScalar sp("sp",dimPressure,1e-2);
    coeffC_ = 0.0*coeff()*neg(p_vap_*xL_ + sp - p_*x_); //no condensation
    coeffV_ = 2.0*betaM_/(2.0-betaM_)*coeff()*pos(p_vap_*xL_ - sp - p_*x_);
    
    m_evap_ = (coeffC_ + coeffV_) * (p_vap_*xL_ - p_*x_);
    
    /*
    m_evap_ = dimensionedScalar("m0",dimMass/dimArea/dimTime,0.1)*pos(area_-sA);
    
    dimensionedScalar totalEvap = fvc::domainIntegrate(m_evap_*area_);
    dimensionedScalar totalArea = fvc::domainIntegrate(area_);
    dimensionedScalar totalLiq = fvc::domainIntegrate(alphaL_);
    
    scalar radius = Foam::sqrt(totalLiq.value()/0.0002/pi)*1000.0;
    
    //volScalarField rhoE = rho_evap(); //m_evap_ * area_
    Foam::Info<< "Net domain evaporation = " << totalEvap.value() << " kg/s" << endl;
    Foam::Info<< "Total area = " << totalArea.value() << " m2" << endl;
    Foam::Info<< "Radius = " << radius << " mm" << endl;
    */
    Foam::Info<< "Min,max evaporation flux for " << vapor_specie_ << " = "
              << Foam::min(m_evap_).value() << ", " 
              << Foam::max(m_evap_).value() << " kg/m2/s" << Foam::endl;
}


// get the explicit and implicit source terms for Yvapor
Pair<tmp<volScalarField> > Foam::evaporationModels::HertzKnudsenPressure::YSuSp() const
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
Pair<tmp<volScalarField> > Foam::evaporationModels::HertzKnudsenPressure::pSuSp() const
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

tmp<volScalarField> Foam::evaporationModels::HertzKnudsenPressure::dPvdT() const
{
    tmp<volScalarField> tdPvdT
    (
        new volScalarField
        (
            IOobject
            (
                "tdPvdT",
                alphaL_.mesh().time().timeName(),
                alphaL_.mesh()
            ),
            alphaL_.mesh(),
            dimensionedScalar("dPdT",dimPressure/dimTemperature,0.0)
        )
    );

    tdPvdT() = p_vap_/T_/T_*
               (
                   Lb_/R_*neg(T_ - Tb_)
                 + (
                       dimensionedScalar("B",dimTemperature,PvCoeffs_[1])
                     + 2*dimensionedScalar("C",dimTemperature*dimTemperature,PvCoeffs_[2])/T_
                   )*pos(T_ - Tb_)
               );

    return tdPvdT;
}

Pair<tmp<volScalarField> > Foam::evaporationModels::HertzKnudsenPressure::TSuSp() const
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
