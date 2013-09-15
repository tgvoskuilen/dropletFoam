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

#include "ThermalDecompReaction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mixturePhaseChangeModels
{
    defineTypeNameAndDebug(ThermalDecompReaction, 0);
    addToRunTimeSelectionTable(mixturePhaseChangeModel, ThermalDecompReaction, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::mixturePhaseChangeModels::ThermalDecompReaction::ThermalDecompReaction
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
    xL_
    (
        IOobject
        (
            "xL_" + name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("xL",dimless,0.0)
    ),
    p_vap_
    (
        IOobject
        (
            "p_vap_" + name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("pv",dimPressure,0.0)
    ),
    Pc_(phaseChangeDict_.lookup("Pc")),
    Tc_(phaseChangeDict_.lookup("Tc")),
    Tb_(phaseChangeDict_.lookup("Tb")),
    Lb_(phaseChangeDict_.lookup("Lb")),
    dH_(phaseChangeDict_.lookup("dH")),
    PvCoeffs_(phaseChangeDict_.lookup("PvCoeffs")),
    betaTD_(phaseChangeDict_.lookup("betaTD")),
    W_(dimensionedScalar("W", dimMass/dimMoles, 0.0)) // kg/kmol
{    
    List<word> prodList = products_.toc();
    List<word> reacList = reactants_.toc();
    
    //TODO: Warn if length of reactants > 1
    
    
    
    //vapor_specie_ = prodList[0];
    liquid_specie_ = reacList[0];
    
    //W_ = alphaV_.subSpecies()[vapor_specie_]->W();
    W_.value() = reacThermo_[liquid_specie_]->W();
    
    Info<< "Liquid decomposition reaction configured for " 
        << liquid_specie_ << endl;
}
    
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::mixturePhaseChangeModels::ThermalDecompReaction::calculate
(
    const volScalarField& phaseChangeZones,
    const volScalarField& area
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
    
    dimensionedScalar sA("sA",dimless/dimLength,SMALL);
    
    // Only decompose above Tb
    mask_ = pos(p_vap_ - p_) * pos(area-sA) * phaseChangeZones;
                        
    //Calculate the mole fraction of liquid specie
    xL_ = alphaL_.x(liquid_specie_);
    
    omega_ = area * betaTD_ * mask_ * xL_ * (p_vap_ - p_);

    Foam::Info<< "Min,max thermal decomposition flux for " << liquid_specie_ << " = "
              << Foam::min(omega_).value() << ", " 
              << Foam::max(omega_).value() << " kmol/m3/s" << Foam::endl;
}

tmp<volScalarField> Foam::mixturePhaseChangeModels::ThermalDecompReaction::mdot
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
    
    if( phaseName == "Liquid" )
    {
        tmdot() -= omega_*W_;
    }
    else if( phaseName == "Vapor" )
    {
        tmdot() += omega_*W_;
    }
    else
    {
        Info<<"WARNING: Invalid phase " << phaseName << endl;
    }
    
    return tmdot;
}

tmp<volScalarField> Foam::mixturePhaseChangeModels::ThermalDecompReaction::Vdot
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
    
    if( phaseName == "Liquid" )
    {
        dimensionedScalar rhoL0 = alphaL_.subSpecies()[liquid_specie_]->rho0();
        tVdot() -= omega_*W_/rhoL0;
    }
    else if( phaseName == "Vapor" )
    {
        //tVdot() += omega_*R_*T_/p_;
        
        scalar sumF = 0.0;
        
        forAllConstIter(HashTable<scalar>, products_, fpI)
        {
            sumF += fpI();
        }
        
        tVdot() += omega_*R_*T_/p_ * sumF;
    }
    else
    {
        Info<<"WARNING: Invalid phase " << phaseName << endl;
    }
    
    return tVdot;
}


// get the explicit and implicit source terms for Yvapor
Pair<tmp<volScalarField> > Foam::mixturePhaseChangeModels::ThermalDecompReaction::YSuSp
(
    const word& specie
) const
{
    Pair<tmp<volScalarField> > tYSuSp
    (
        tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "tYSu",
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar("YSu", dimDensity/dimTime, 0.0)
            )
        ),
        tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "tYSp",
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar("YSp", dimDensity/dimTime, 0.0)
            )
        )
    );
    
    
    /*if( specie == vapor_specie_ )
    {
        //S_Yv = explicit - implicit*Yv
        // This casts the evaporation in the form C*(Ysat - Y)
         
        tmp<volScalarField> xByY = alphaV_.xByY(vapor_specie_);
        tYSuSp.first() = (coeffC_+coeffV_)*p_vap_*xL_;
        tYSuSp.second() = (coeffC_+coeffV_)*xByY*p_;
    }*/
    
    if( specie == liquid_specie_ )
    {
        //S_YL = explicit
        //tYSuSp.first() = -(coeffC_+coeffV_)*(p_vap_*xL_ - x_*p_);
        tYSuSp.first() = -omega_*W_;
    }
    else if( products_.found(specie) )
    {
        //scalar Wp = prodThermo_[specie]->W()/W_.value() * products_[specie];
        
        dimensionedScalar Wp("Wp",dimMass/dimMoles,prodThermo_[specie]->W());
        tYSuSp.first() = omega_ * Wp * products_[specie];
        
        /*tmp<volScalarField> xByY = alphaV_.xByY(specie);
        tYSuSp.first() = (coeffC_+coeffV_)*p_vap_*xL_ * Wp;
        tYSuSp.second() = (coeffC_+coeffV_)*xByY*p_ * Wp;*/
    }

    return tYSuSp;
}




Pair<tmp<volScalarField> > Foam::mixturePhaseChangeModels::ThermalDecompReaction::TSuSp() const
{
    Pair<tmp<volScalarField> > tTSuSp
    (
        tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "tTSu",
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar("TSu", dimPower/dimVolume, 0.0)
            )
        ),
        tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "tTSp",
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar("TSp", dimPower/dimVolume/dimTemperature, 0.0)
            )
        )
    );

    //Sh = -omega_*L()
    //
    //Sh = explicit - implicit * T
    //
    //  implicit = -dSh/dT = (dm/dT*L() + m_evap*dL/dT)*area_
    //    dm/dT = dcoeff/dT*(pv*xL-p*x) + coeff*xL*dpv/dT
    //    dcoeff/dT = -coeff/(2*T)
    //    dpv/dT = pv * L / (R*T^2)

    /*tmp<volScalarField> tL = L();
    tmp<volScalarField> dodT = -omega_/(2*T_) + (coeffC_+coeffV_)/W_*xL_*dPvdT();
    tmp<volScalarField> Sp = dodT*tL() + omega_*dLdT();
    tmp<volScalarField> Su = -omega_*tL() + Sp()*T_;
    
    return Pair<tmp<volScalarField> >
    (
        Su,
        Sp
    );*/
    
    tTSuSp.first() = omega_ * dH_;
    
    return tTSuSp;

}




// ************************************************************************* //
