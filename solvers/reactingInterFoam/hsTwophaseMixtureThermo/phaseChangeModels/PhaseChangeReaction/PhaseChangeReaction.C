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
    tmp<volScalarField> xByY = alphaV_.xByY(vapor_specie_);
    
    //S_Yv = explicit - implicit*Yv
    // This casts the evaporation in the form C*(Ysat - Y)
        
    return Pair<tmp<volScalarField> >
    (
        area_*(coeffC_+coeffV_)*p_vap_*xL_,
        area_*(coeffC_+coeffV_)*xByY*p_
    );
}



Pair<tmp<volScalarField> > Foam::phaseChangeModels::PhaseChangeReaction::TSuSp() const
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
    
    
    return tVdot;
}
// ************************************************************************* //
