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
    rhoPhiAlpha_
    (
        IOobject
        (
            "rhoPhiAlpha"+name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rhoPhiAlpha"+name, dimMass/dimTime, 0.0)
    ),
    rhoAlpha_
    (
        IOobject
        (
            "rhoAlpha_"+name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rhoAlpha_"+name, dimDensity, 1.0)
    ),
    otherPhase_(NULL),
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
    ),
    cellMask_
    (
        IOobject
        (
            "cellMask_"+name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("cellMask_"+name, dimless, 0.0)
    ),
    faceMask_
    (
        IOobject
        (
            "faceMask_"+name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("faceMask_"+name, dimless, 0.0)
    )
{  
    rhoAlpha_.oldTime();
    if( phaseDict_.found("SchmidtNo") )
    {
        Sc_ = readScalar(phaseDict_.lookup("SchmidtNo"));
    }
    
    Info<< "Created phase " << name << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phase> Foam::phase::clone() const
{
    notImplemented("phase::clone() const");
    return autoPtr<phase>(NULL);
}

void Foam::phase::setSpecies( const volScalarField& rhoTotal )
{
    forAllIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {   
        specieI().Yp() = specieI().Y() / (Yp() + SMALL);
        specieI().Yp().oldTime();
    }
}

// Only applicable for liquid phase. Vapor phase will return 0 here.
// TODO: Have vapor specie do its sutherland calc here too
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
            tmu() += specieI().nuModel().nu() * specieI().Yp() * rho(p,T);
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


void Foam::phase::correct(const volScalarField& p, const volScalarField& T)
{
    rhoAlpha_ = (*this) * rho(p,T);
    
    Info<< "Min,max rhoAlpha"<<name_
        <<" = " << Foam::min(rhoAlpha_).value() << ", " 
        << Foam::max(rhoAlpha_).value() << endl;
        
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
        forAllConstIter
        (
            PtrDictionary<subSpecie>, 
            otherPhase_->subSpecies(), 
            specieLI
        )
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

    dimensionedScalar dA = Foam::pow(Foam::min(mesh().V()),2.0/3.0)/30.0;
    
    forAllConstIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {
        if( specieI().sigma0().value() > SMALL )
        {
            volScalarField Ys(specieI().Y());
            
            for(label i = 0; i < 5; i++)
            {
                Ys += dA * fvc::laplacian(Ys);
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


Foam::tmp<volScalarField> Foam::phase::Ypp() const
{
    tmp<volScalarField> tYpp
    (
        new volScalarField
        (
            IOobject
            (
                "tYpp"+name_,
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("Ypp",dimless,0.0)
        )
    );
    
    forAllConstIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {
        tYpp() += specieI().Yp();
    }
    
    return tYpp;
}


// Volumetric sink term from evaporation = m_evap/rhoL (only applicable to
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
            tSv() += specieI().evapModel().m_evap() / specieI().rho0();
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
        return psi(T) * p;
    }
    else
    {
        dimensionedScalar rhoBase("rhoBase",dimDensity,1000);
        
        tmp<volScalarField> Yp_ = Yp();
        tmp<volScalarField> Yvoid = 0.0001*neg(Yp_()-0.05);
        tmp<volScalarField> den = Yvoid()/rhoBase;
        
        forAllConstIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
        {
            den() += specieI().Y() / specieI().rho0();
        }
        
        return (Yvoid + Yp_)/den;
    }
}


// Construct the compressibility field of the phase
// for Vapors, psi = W/(R*T)
// for liquids, psi = 0
Foam::tmp<volScalarField> Foam::phase::psi
(
    const volScalarField& T
) const
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

    if (name_ == "Vapor")
    {
        dimensionedScalar Ru("Ru",dimensionSet(1, 2, -2, -1, -1),8314);
        tpsi() = W() / (Ru * T);
    }

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


// Construct the molar mass field.
Foam::tmp<volScalarField> Foam::phase::W() const
{
//    return Yp()/Np();
    dimensionedScalar ws("ws",dimMass/dimMoles,SMALL);
    tmp<volScalarField> Wother = otherPhase_->Yp() / otherPhase_->Np();
    
    tmp<volScalarField> Yp_ = Yp();
    tmp<volScalarField> Yvoid = 0.0001*neg(Yp_()-0.05);
    tmp<volScalarField> den = Yvoid()/(Wother + ws);
    
    forAllConstIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {
        den() += specieI().Y() / specieI().W();
    }
    
    return (Yvoid + Yp_)/den;
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


void Foam::phase::calculateDs
(
    const volScalarField& muEff
)
{
    D_ = Sc_ * fvc::interpolate(muEff) * faceMask_;
}


/*tmp<surfaceScalarField> Foam::phase::calculateDs
(
    const compressible::turbulenceModel& turb
)
{
    Info<< "Calculating diffusion coefficients for " << name_ << endl;
    
    tmp<surfaceScalarField> tUc
    (
        new surfaceScalarField
        (
            IOobject
            (
                "tUc"+name_,
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("Uc",dimMass/dimTime,0.0)
        )
    );

    D_ = Sc_ * fvc::interpolate(turb.muEff()) * faceMask_;
    
    forAllIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {
        const volScalarField& Yi = specieI().Yp();
        tUc() += D_ * fvc::snGrad(Yi) * mesh().magSf();
    }

    return tUc;
}*/


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

// Returns a sharpened version of the phase volume fraction
Foam::tmp<Foam::volScalarField> Foam::phase::sharp
(
    scalar tol
) const
{
    tmp<volScalarField> ts
    (
        new volScalarField
        (
            IOobject
            (
                "tAlphaSharp"+name_,
                mesh().time().timeName(),
                mesh()
            ),
            *this
        )
    );

    volScalarField& s = ts();

    //Sharpen the alpha field
    scalar Cpc = 2.0*tol;
    s = (Foam::min(Foam::max(s, 0.5*Cpc),1.0-0.5*Cpc) - 0.5*Cpc)/(1.0-Cpc);

    return ts;
}

        
// Solve the mass fraction governing equation (YEqn) for each subspecie of
// this phase.
scalar Foam::phase::solveSubSpecies
(
    const volScalarField& p,
    const volScalarField& T
)
{   
    word divScheme("div(rho*phi*alpha,Yi)");

    tmp<volScalarField> DbyRho = fvc::average(D_) / rho(p,T);
    scalar MaxFo = mesh().time().deltaTValue() * Foam::max
    (
        DbyRho().internalField() / pow(mesh().V(),2.0/3.0)
    );
        
    Info<< "Max D = " << Foam::max(D_).value() << endl;
    Info<< "Max Fo = " << MaxFo << endl;
        
    Info<< "Solving "<<subSpecies_.size()
        <<" subspecie(s) for phase: " 
        << name_ << endl;
        
    //solve ddt(rhoAlpha) + div(rhoPhiAlpha) == 0
    scalarField& rhoAlphaIf = rhoAlpha_;
    const scalarField& rhoAlpha0 = rhoAlpha_.oldTime();
    const scalar deltaT = mesh().time().deltaTValue();

    rhoAlphaIf = 0.0;
    surfaceScalarField rhoPhiAlphaMasked = rhoPhiAlpha_*faceMask_;
    fvc::surfaceIntegrate(rhoAlphaIf, rhoPhiAlphaMasked);
    
    //TODO: Need source term in here for evaporation
    forAll(rhoAlphaIf, cellI)
    {
        rhoAlphaIf[cellI] = rhoAlpha0[cellI] - rhoAlphaIf[cellI]*deltaT;
    }
    rhoAlpha_.correctBoundaryConditions();
    
    volScalarField Sp = (1.0 - cellMask_)*rho(p,T)/mesh().time().deltaT();
    rhoAlpha_ *= cellMask_;
    
   
    // Loop through phase's subspecies
    forAllIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {
        volScalarField& Yi = specieI().Yp();
        
        fvScalarMatrix YiEqn
        (
            fvm::ddt(rhoAlpha_, Yi)
          + fvm::div(rhoPhiAlphaMasked, Yi, divScheme)
          - fvm::laplacian(D_, Yi)
         ==
            combustionPtr_->R(Yi)
          + Su_Yi_evap( specieI() )
          - fvm::Sp(Sp, Yi)
        );

        YiEqn.relax();
        YiEqn.solve(mesh().solver("Yi"));
        
        Yi.max(0.0);
        Yi.min(1.0);
        
        Info<< Yi.name() << " min,max,avg = " 
            << Foam::min(Yi).value() <<  ", " << Foam::max(Yi).value() 
            << ", " << Yi.weightedAverage(mesh().V()).value() << endl;
            
    }
        
    return MaxFo;
}

void Foam::phase::setPhaseMasks(scalar maskTol)
{
    const volScalarField& alpha = *this;
    cellMask_ = pos(alpha - maskTol);
    faceMask_ = pos(fvc::interpolate(cellMask_) - 0.95);
}

void Foam::phase::updateGlobalYs(const volScalarField& rhoAlphaOther)
{
    forAllIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {   
        specieI().Y() = specieI().Yp()*rhoAlpha_/(rhoAlpha_+rhoAlphaOther);
        specieI().Y().correctBoundaryConditions();
    }
}


// ************************************************************************* //
