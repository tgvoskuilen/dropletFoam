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

#include "phase.H"
#include "subSpecie.H"
#include "evaporationModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phase::phase
(
    const word& name,
    const dictionary& phaseDict,
    const fvMesh& mesh,
    PtrList<volScalarField>& species,
    const PtrList<gasThermoPhysics>& speciesData
)
:
    volScalarField
    (
        IOobject
        (
            "alpha" + name,
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    name_(name),
    phaseDict_(phaseDict),
    rhoPhi_
    (
        IOobject
        (
            "rhoPhi"+name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rhoPhi"+name, dimMass/dimTime, 0.0)
    ),
    rho_
    (
        IOobject
        (
            "rho"+name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rho"+name, dimDensity, 1.2)
    ),
    T_
    (
        IOobject
        (
            "T"+name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh.lookupObject<volScalarField>("T")
    ),
    mu_
    (
        IOobject
        (
            "mu"+name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("mu",dimMass/dimLength/dimTime,0.0)
        //mesh.lookupObject<volScalarField>("mu")
    ),
    combustionPtr_(NULL),
    species_(species),
    speciesData_(speciesData),
    subSpecies_
    (
        phaseDict_.lookup("subspecies"),
        subSpecie::iNew(mesh, species, speciesData)
    ),
    Sc_(1.0),
    D_
    (
        IOobject
        (
            "D_"+name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("D_"+name, dimDensity*dimArea/dimTime, 0.0)
    )
{    
    if( name_ == "Liquid" )
    {
        mu_ = dimensionedScalar("mu",dimMass/dimLength/dimTime,1e-3);
    }
    else
    {
        mu_ = dimensionedScalar("mu",dimMass/dimLength/dimTime,1.8e-5);
    }

    rho_ = rho(T_,T_); //arguments not used for now
    rho_.oldTime();
    
    T_.oldTime();
    
    if( phaseDict_.found("SchmidtNo") )
    {
        Sc_ = readScalar(phaseDict_.lookup("SchmidtNo"));
    }
        
    Foam::Info<< "Created phase " << name << Foam::endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phase> Foam::phase::clone() const
{
    notImplemented("phase::clone() const");
    return autoPtr<phase>(NULL);
}

/*
Foam::tmp<Foam::volScalarField> Foam::phase::U_Sp
(
    const PLICInterface& i
)
{
    if( name_ == "Vapor" )
    {
        dimensionedScalar rhorDT
        (
            "rhorDT",
            dimDensity/dimTime,
            1.0/mesh().time().deltaTValue()
        );
        U_Sp_ = (i.smallGasCells() + i.noGasCells())*rhorDT;
        return (i.smallGasCells() + i.noGasCells())*rhorDT;
    }
    else
    {
        dimensionedScalar rhorDT
        (
            "rhorDT",
            dimDensity/dimTime,
            1000./mesh().time().deltaTValue()
        );
        U_Sp_ = (i.smallLiquidCells() + i.noLiquidCells())*rhorDT;
        return (i.smallLiquidCells() + i.noLiquidCells())*rhorDT;
    }
}*/
/*
Foam::tmp<Foam::volVectorField> Foam::phase::U_Su
(
    const PLICInterface& i,
    const volVectorField& U_opp,
    const volScalarField& mu_opp
)
{
    if( name_ == "Vapor" )
    {
        dimensionedScalar rhorDT
        (
            "rhorDT",
            dimDensity/dimTime,
            1.0/mesh().time().deltaTValue()
        );
        
        U_Su_ = (i.smallGasCells()*i.scAverage<vector>(name_,U_)
             + i.noGasCells()*UFarField(i)) * rhorDT
             + i.gasCells()*i.shearVec(name_, U_opp, U_, mu_opp, mu_);
        
        return (i.smallGasCells()*i.scAverage<vector>(name_,U_)
             + i.noGasCells()*UFarField(i)) * rhorDT
             + i.gasCells()*i.shearVec(name_, U_opp, U_, mu_opp, mu_);
    }
    else
    {
        dimensionedScalar rhorDT
        (
            "rhorDT",
            dimDensity/dimTime,
            1000.0/mesh().time().deltaTValue()
        );
        U_Su_ = (i.noLiquidCells() * UFarField(i)
             + i.smallLiquidCells()*i.scAverage<vector>(name_,U_)) * rhorDT
             + i.liquidCells()*i.shearVec(name_, U_, U_opp, mu_, mu_opp);
        
        return (i.noLiquidCells() * UFarField(i)
             + i.smallLiquidCells()*i.scAverage<vector>(name_,U_)) * rhorDT
             + i.liquidCells()*i.shearVec(name_, U_, U_opp, mu_, mu_opp);
    }
}*/



void Foam::phase::setRhoPhi
(
    const surfaceScalarField& phiAlpha,
    scalar fraction
)
{
    //TODO: Correct interpolation near interface
    if( fraction > SMALL )
    {
        rhoPhi_ += fraction * phiAlpha * fvc::interpolate(rho_);
    }
    else
    {
        rhoPhi_ = phiAlpha * fvc::interpolate(rho_);
    }
}


Foam::tmp<Foam::volScalarField> Foam::phase::alphaCorr
(
    const PLICInterface& i
) const
{
    const volScalarField& alpha = *this;
    return alpha * pos(alpha - i.alphaMin());
}



label Foam::phase::numSmallAndEmptyCells
(
    const PLICInterface& i
) const
{
    tmp<volScalarField> mask;
    if( name_ == "Vapor" )
    {
        mask = i.smallGasCells() + i.noGasCells();
    }
    else
    {
        mask = i.smallLiquidCells() + i.noLiquidCells();
    }
    
    return label(Foam::sum(mask()).value());
}


/*
void Foam::phase::calcFixedRhoValues
(
    labelList& cells,
    scalarList& values,
    const PLICInterface& i
) const
{
    tmp<volScalarField> rhovalue;
    tmp<volScalarField> mask;
    if( name_ == "Vapor" )
    {
        mask = i.smallGasCells() + i.noGasCells();
        rhovalue = i.noGasCells()*rhoFarField(i)
                 + i.smallGasCells()*i.scAverage<scalar>(name_,rho_);
    }
    else
    {
        mask = i.smallLiquidCells() + i.noLiquidCells();
        rhovalue = i.noLiquidCells()*rhoFarField(i)
                 + i.smallLiquidCells()*i.scAverage<scalar>(name_,rho_);
    }
        
    label idx = 0;
    forAll(mask(), cellI)
    {
        if( mask()[cellI] > 0.5 )
        {
            cells[idx] = cellI;
            values[idx] = rhovalue()[cellI];
            ++idx;
        }
    }
}
*/


void Foam::phase::updateRho
(
    const PLICInterface& interface
)
{
    const volScalarField& alpha = *this;
    Info<< "Max rhoErr" << name_ << " = " 
        << Foam::max(Foam::mag(fvc::ddt(alpha,rho_) + fvc::div(rhoPhi_)))
        << endl;

/*
    volScalarField rhoTmp = rho_;
    
    Foam::solve( fvm::ddt(alpha,rhoTmp) + fvc::div(rhoPhi_), mesh().solverDict("rho") );
    
    Info<< "Min,max rhoTmp" << name_ << " = " 
        << Foam::min(rhoTmp).value() << ", " << Foam::max(rhoTmp).value() <<
        << endl;
*/
}

/*
Foam::tmp<Foam::volScalarField> Foam::phase::rho_Sp
(
    const PLICInterface& i
) const
{
    dimensionedScalar rDT = 1.0 / mesh().time().deltaT();
    
    if( name_ == "Vapor" )
    {
        return (i.smallGasCells() + i.noGasCells())*rDT;
    }
    else
    {
        return (i.smallLiquidCells() + i.noLiquidCells())*rDT;
    }
}


Foam::dimensionedScalar Foam::phase::rhoFarField(const PLICInterface& i) const
{
    return rho_.weightedAverage(i.area());
}

Foam::dimensionedVector Foam::phase::UFarField(const PLICInterface& i) const
{
    return U_.weightedAverage(i.area());
}

Foam::dimensionedScalar Foam::phase::pFarField(const PLICInterface& i) const
{
    return p_.weightedAverage(i.area());
}



Foam::tmp<Foam::volScalarField> Foam::phase::rho_Su
(
    const PLICInterface& i
) const
{
    dimensionedScalar rDT = 1.0 / mesh().time().deltaT();
    
    if( name_ == "Vapor" )
    {
        return (i.noGasCells() * rhoFarField(i)
             + i.smallGasCells()*i.scAverage<scalar>(name_,rho_))*rDT; //Add m_evap
    }
    else
    {
        return (i.noLiquidCells() * rhoFarField(i)
             + i.smallLiquidCells()*i.scAverage<scalar>(name_,rho_))*rDT; //Add m_evap
    }
}




Foam::tmp<Foam::volScalarField> Foam::phase::p_Sp
(
    const PLICInterface& i,
    const volScalarField& rhorAU
)
{
    
    if( name_ == "Vapor" )
    {
        dimensionedScalar psirDT
        (
            "psirDT",
            dimTime/dimArea,
            1e-5/mesh().time().deltaTValue()
        );
        p_Sp_ = (i.smallGasCells() + i.noGasCells())*psirDT;
        return (i.smallGasCells() + i.noGasCells())*psirDT;
    }
    else
    {
        dimensionedScalar psirDT
        (
            "psirDT",
            dimTime/dimArea,
            1.0/mesh().time().deltaTValue()
        );
        
        p_Sp_ = i.noLiquidCells()*psirDT
             + i.smallLiquidCells()*psirDT
             + i.liquidCells()*rhorAU*i.AbydeltaLV();
             
        return i.noLiquidCells()*psirDT
             + i.smallLiquidCells()*psirDT
             + i.liquidCells()*rhorAU*i.AbydeltaLV();
    }
}


Foam::tmp<Foam::volScalarField> Foam::phase::p_Su
(
    const PLICInterface& i,
    const volScalarField& p_opp,
    const volScalarField& rhorAU
)
{
    if( name_ == "Vapor" )
    {
        dimensionedScalar psirDT
        (
            "psirDT",
            dimTime/dimArea,
            1e-5/mesh().time().deltaTValue()
        );
        dimensionedScalar p0("p0",dimPressure,1e5);
        
        p_Su_ = (i.noGasCells()*pFarField(i)
             + i.smallGasCells()*i.scAverage<scalar>(name_,p_))*psirDT; //Add m_evap
             
        return (i.noGasCells()*pFarField(i)
             + i.smallGasCells()*i.scAverage<scalar>(name_,p_))*psirDT; //Add m_evap
    }
    else
    {
        //TODO: TEMPORARY CONSTANT
        dimensionedScalar sigmaT("sigma",dimPressure*dimLength,0.07);
        
        dimensionedScalar psirDT
        (
            "psirDT",
            dimTime/dimArea,
            1.0/mesh().time().deltaTValue()
        );
        
        p_Su_ = i.noLiquidCells()*pFarField(i)*psirDT
             + i.smallLiquidCells()*(p_opp + sigmaT*i.kappa())*psirDT
             + i.liquidCells()*rhorAU*i.AbydeltaLV()*(p_opp + sigmaT*i.kappa()); //Add m_evap
             
        return i.noLiquidCells()*pFarField(i)*psirDT
             + i.smallLiquidCells()*(p_opp + sigmaT*i.kappa())*psirDT
             + i.liquidCells()*rhorAU*i.AbydeltaLV()*(p_opp + sigmaT*i.kappa()); //Add m_evap
    }
}
*/

// Only applicable for liquid phase. Vapor phase will return 0 here.
Foam::tmp<Foam::volScalarField> Foam::phase::mu
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    tmp<volScalarField> tmu
    (
        new volScalarField
        (
            IOobject
            (
                "tmu",
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("tmu", dimArea*dimDensity/dimTime, 0.0)
        )
    );
    
    forAllConstIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {   
        if( specieI().hasNuModel() )
        {
            tmu() += specieI().nuModel().nu() * specieI().Y() * rho(p,T);
        }
    }
    
    return tmu;
}



Foam::volScalarField& Foam::phase::Y(const word& specie)
{
    return subSpecies_[specie]->Y();
}

const Foam::volScalarField& Foam::phase::Y(const word& specie) const
{
    return subSpecies_[specie]->Y();
}


void Foam::phase::setCombustionPtr
(
    combustionModels::rhoChemistryCombustionModel* combustion
)
{
    combustionPtr_ = combustion;
}


void Foam::phase::correct()
{
    //rho_ = rho(p_,T_);

    forAllIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {   
        specieI().correct();
    }
}


bool Foam::phase::read(const dictionary& phaseDict)
{
    phaseDict_ = phaseDict;

    return true;
}

// Calculates the total net volumetric source due to evaporation. Only relevant
// for the liquid phase, and called by S_evap() in thermo for in pEqn
Foam::tmp<Foam::volScalarField> Foam::phase::S_evap
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    tmp<volScalarField> tS_evap
    (
        new volScalarField
        (
            IOobject
            (
                "tS_evap",
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("S_evap", dimless/dimTime, 0.0)
        )
    );
    
    dimensionedScalar Ru("Ru", dimEnergy/dimMoles/dimTemperature, 8314);
    
    forAllConstIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {
        dimensionedScalar rho0 = specieI().rho0();
        
        if( specieI().hasEvaporation() && rho0.value() > SMALL )
        {
            dimensionedScalar Wv = specieI().W();
            
            tS_evap() += specieI().evapModel().m_evap()*(Ru*T/(p*Wv) - 1/rho0);
        }
    }
    
    return tS_evap;
}



// Mass source term due to evaporation for a given subspecie. This is a source
// term in the Yi equation for each subspecie
Foam::tmp<Foam::volScalarField> Foam::phase::Su_Yi_evap
(
    const phase& alphaL,
    const subSpecie& specieI
) const
{
    tmp<volScalarField> tSu_evap
    (
        new volScalarField
        (
            IOobject
            (
                "tS_Yiu_evap",
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("S_Yiu_evap", dimDensity/dimTime, 0.0)
        )
    );
    
    if (name_ == "Vapor")
    {
        forAllConstIter(PtrDictionary<subSpecie>, alphaL.subSpecies(), specieLI)
        {
            if (specieLI().hasEvaporation())
            {
                if (specieLI().evapModel().vaporName() == specieI.Y().name())
                {
                    tSu_evap() += specieLI().evapModel().m_evap();
                    break;
                }
            }
        }
    }
    else
    {
        if (specieI.hasEvaporation())
        {
            tSu_evap() -= specieI.evapModel().m_evap();
        }
    }
    
   
    return tSu_evap;
}



Foam::tmp<volScalarField> Foam::phase::sigma
(
    const volScalarField& T,
    const volScalarField& kappaMask
) const
{
    tmp<volScalarField> tsigma
    (
        new volScalarField
        (
            IOobject
            (
                "tsigma"+name_,
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("sigma", dimensionSet(1, 0, -2, 0, 0), 0.0)
        )
    );
    
    tmp<volScalarField> tn
    (
        new volScalarField
        (
            IOobject
            (
                "tn"+name_,
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("n", dimless, 0.0)
        )
    );
    
    dimensionedScalar dTau("dTau",dimArea,1e-12);
    
    forAllConstIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {
        if( specieI().sigma0().value() > SMALL )
        {
            volScalarField Ys(specieI().Y());
            
            for(label i = 0; i < 3; i++)
            {
                Ys += dTau * fvc::laplacian(Ys);
            }
            
            Ys = pos(Foam::mag(Ys) - SMALL) * kappaMask;
            tn() += Ys;
            tsigma() += Ys*specieI().sigma(T);
        }
    }
    
    tn() += neg(tn() - SMALL);

    return tsigma/tn;
}



// calculate the net heat source/sink for evaporation for this phase. This is
// called by thermo when assembling the total sources for the TEqn. It
// is only called on the liquid phase.
Foam::tmp<volScalarField> Foam::phase::Sh_evap() const
{
    tmp<volScalarField> tSh
    (
        new volScalarField
        (
            IOobject
            (
                "tSh"+name_,
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("tSh"+name_, dimPower/dimVolume, 0.0)
        )
    );
    
    // Get evaporation enthalpy source
    forAllConstIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {
        if(specieI().hasEvaporation())
        {
            tSh() += specieI().evapModel().Sh();
        }
    }
    
    return tSh;    
}


// Get the mass fraction sum of this phase
Foam::tmp<volScalarField> Foam::phase::Yp() const
{
    tmp<volScalarField> tYp
    (
        new volScalarField
        (
            IOobject
            (
                "tYp"+name_,
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("Yp",dimless,0.0)
        )
    );
    
    forAllConstIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {
        tYp() += specieI().Y();
    }
    
    return tYp;
}


// Volumetric source term from evaporation = -m_evap/rhoL (only applicable to
//  liquid phase, used in the alpha equation in the source term)
Foam::tmp<volScalarField> Foam::phase::Sv_evap() const
{
    tmp<volScalarField> tSv
    (
        new volScalarField
        (
            IOobject
            (
                "tSv"+name_,
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("tSv"+name_, dimless/dimTime, 0.0)
        )
    );

    forAllConstIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {
        if( specieI().hasEvaporation() )
        {
            tSv() -= specieI().evapModel().m_evap() / specieI().rho0();
        }
    }

    return tSv;
}


// Construct the density field of the phase. For the vapor phase will never
// return a zero density (rho = p * W / (R * T)) unless W is zero
// the solid phase should not return zeros, but may return very large #s?
Foam::tmp<volScalarField> Foam::phase::rho
(
    const volScalarField& p, 
    const volScalarField& T
) const
{
    tmp<volScalarField> trho
    (
        new volScalarField
        (
            IOobject
            (
                "trho"+name_,
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("trho"+name_, dimDensity, 0.0)
        )
    );
    
    if (name_ == "Vapor")
    {
        trho() = dimensionedScalar("rhoV",dimDensity,1.2);
        return trho;
        
        //return psi() * p;
    }
    else
    {
        trho() = dimensionedScalar("rhoL",dimDensity,1000.0);
        return trho;
    /*
        tmp<volScalarField> den
        (
            new volScalarField
            (
                IOobject
                (
                    "den",
                    mesh().time().timeName(),
                    mesh()
                ),
                mesh(),
                dimensionedScalar("den", dimless/dimDensity, SMALL)
            )
        );
        
        forAllConstIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
        {
            den() += specieI().Y() / specieI().rho0();
        }
        
        return Yp() / den;
        
        */
    }
}


// Construct the compressibility field of the phase
// for Vapors, psi = W/(R*T)
// for liquids, psi = 0
Foam::tmp<volScalarField> Foam::phase::psi() const
{
    tmp<volScalarField> tpsi
    (
        new volScalarField
        (
            IOobject
            (
                "tpsi"+name_,
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar
            (
                "tpsi"+name_, 
                dimTime*dimTime/dimLength/dimLength, 
                0.0
            )
        )
    );

    /*if (name_ == "Vapor")
    {
        tpsi() = W()/(dimensionedScalar("R", dimensionSet(1, 2, -2, -1, -1), 8314) * T_);
    }*/

    return tpsi;
}


Foam::tmp<volScalarField> Foam::phase::Np() const
{
    //Np = sum(Yi/Wi)
    tmp<volScalarField> tNp
    (
        new volScalarField
        (
            IOobject
            (
                "tNp"+name_,
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("tNp"+name_, dimensionSet(-1,0,0,0,1), SMALL)
        )
    );
    
    forAllConstIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {
        tNp() += specieI().Y() / specieI().W();
    }
    
    return tNp;
}


// Construct the molar mass field. Should be zero if Yp == 0. If Yp == VSMALL
// though then W -> VSMALL
Foam::tmp<volScalarField> Foam::phase::W() const
{
    return Yp()/Np();
}



// Calculate the phase thermal conductivity using the method from Harvazinski
// This is only used for the liquid phase. The gas phase conducitivty is 
// evaluated using the selected transport model (Sutherland)
Foam::tmp<Foam::volScalarField> Foam::phase::k() const
{
    tmp<volScalarField> tk1
    (
        new volScalarField
        (
            IOobject
            (
                "tk1"+name_,
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("k1",dimPower/dimLength/dimTemperature,0)
        )
    );
    
    tmp<volScalarField> tk2
    (
        new volScalarField
        (
            IOobject
            (
                "tk2"+name_,
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("k2",dimLength*dimTemperature/dimPower,SMALL)
        )
    );
    
    forAllConstIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {
        // xi = Y / (W * Np)
        tk1() += specieI().Y() / (specieI().W() * Np()) * specieI().kappa();
        tk2() += specieI().Y() / (specieI().W() * Np()) / specieI().kappa();
    }
    
    tk1() += pos(Yp() - 1e-3) / tk2;

    return 0.5*tk1;
}



tmp<surfaceVectorField> Foam::phase::calculateDs
(
    bool allowDiffusion,
    const compressible::turbulenceModel& turb,
    const surfaceScalarField& alphaf
)
{
    Info<< "Calculating diffusion coefficients for " << name_ << endl;
    
    tmp<surfaceVectorField> tUc
    (
        new surfaceVectorField
        (
            IOobject
            (
                "tUc"+name_,
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedVector("Uc",dimVelocity*dimDensity,vector::zero)
        )
    );
        

    //const phase& alpha = *this;
    tmp<volScalarField> tYp = Yp();
    
    D_ = Sc_ * fvc::interpolate(turb.muEff()) * alphaf;
    
    if( allowDiffusion )
    {
        forAllIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
        {
            const volScalarField& Yi = specieI().Y();
            //tUc() += D_ * fvc::interpolate( fvc::grad(Yi / (tYp()+SMALL)) );
            tUc() += D_ * fvc::interpolate( fvc::grad(Yi) / (tYp()+SMALL) );
        }
    }
    else
    {
        D_ *= 0.0;
    }

    return tUc;
}


// Calculate phase specific heat using subspecie Cv model, which calls the
// underlying thermo model for each specie (Janaf table)
// Where the phase is trace (Yp < 1e-4) the Cv is just the average of the
// phase's subspecies, otherwise it is weighted by mass fraction
Foam::tmp<Foam::volScalarField> Foam::phase::Cv
(
    const volScalarField& T
) const
{
    tmp<volScalarField> tCv
    (
        new volScalarField
        (
            IOobject
            (
                "tCv"+name_,
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("Cv",dimEnergy/dimMass/dimTemperature,SMALL)
        )
    );
    
    label ns = subSpecies_.size();
    
    tmp<volScalarField> tYp = Yp();
    tmp<volScalarField> mask1 = pos(tYp() - 1e-4);
    tmp<volScalarField> mask2 = (1 - mask1()) / scalar(ns);
    mask1() /= (tYp + SMALL);
    
    forAllConstIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {
        tCv() += ( mask1()*specieI().Y() + mask2() )*specieI().Cv(T);
    }
    
    tCv().correctBoundaryConditions();
    
    return tCv;
}

        
// Solve the mass fraction governing equation (YEqn) for each subspecie of
// this phase.
scalar Foam::phase::solveSubSpecies
(
    const volScalarField& rhoTotal,
    const surfaceScalarField& rhoPhi,
    const volScalarField& p,
    const volScalarField& T,
    const phase& alphaLiquid,
    const surfaceVectorField& uc,
    const tmp<fv::convectionScheme<scalar> >& mvConvection,
    const surfaceScalarField& alphaf
)
{   
    word divScheme("div(rho*phi*alpha,Yi)");
    word divSchemeCorr("div(uc*Yi)");
    scalar MaxFo = Foam::max(fvc::average(D_)().internalField()/rhoTotal.internalField()/Foam::pow(mesh().V(),2.0/3.0))*mesh().time().deltaTValue();
    
    Info<< "Max D = " << Foam::max(D_).value() << endl;
    Info<< "Max Fo = " << MaxFo << endl;
    
    const phase& alpha = *this;
        
    Info<< "Solving " << subSpecies_.size()
        <<" subspecie(s) for phase: " << name_ << endl;
        
    // Save Yp so it is constant throughout subspecie solving
    tmp<volScalarField> Yp0 = Yp()+SMALL;
    tmp<surfaceScalarField> Yp0f = fvc::interpolate(Yp0());
    
    // Loop through phase's subspecies
    forAllIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {
        volScalarField& Yi = specieI().Y();
        
        fvScalarMatrix YiEqn
        (
            fvm::ddt(rhoTotal, Yi)
        //+ mvConvection->fvmDiv(rhoPhiAlpha_, Yi)
        //  + fvm::div(rhoPhiAlpha_*alphaf, Yi, divScheme)
          + mvConvection->fvmDiv(rhoPhi, Yi)
         ==
          //  fvm::laplacian(D_/Yp0f(), Yi)
          //- fvm::div((uc/Yp0f()) & mesh().Sf(), Yi, divSchemeCorr)
          //- fvc::laplacian(D_*fvc::interpolate(Yi)/Yp0f()/Yp0f(), Yp0f())
           alpha*combustionPtr_->R(Yi)
          + Su_Yi_evap( alphaLiquid, specieI() )
        );

        YiEqn.relax();
        YiEqn.solve(mesh().solver("Yi"));
        
        Yi.max(0.0);
        Yi.min(1.0);
        
        Info<< "  Pre-coerce: " << Yi.name() << " min,max,avg = " 
            << Foam::min(Yi).value() <<  ", " << Foam::max(Yi).value() 
            << ", " << Yi.weightedAverage(mesh().V()).value() << endl;
    }

    return MaxFo;
}



// ************************************************************************* //
