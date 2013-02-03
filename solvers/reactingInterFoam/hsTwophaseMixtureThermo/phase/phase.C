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
#include "diffusionModel.H"

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
    nuModel_
    (
        viscosityModel::New
        (
            "nu" + name, 
            phaseDict_, 
            this->db().lookupObject<volVectorField>("U"), 
            this->db().lookupObject<surfaceScalarField>("phi")
        )
    ),
    combustionPtr_(NULL),
    species_(species),
    speciesData_(speciesData),
    subSpecies_
    (
        phaseDict_.lookup("subspecies"),
        subSpecie::iNew(mesh, species, speciesData)
    )
{
    Foam::Info<< "Created phase " << name << Foam::endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phase> Foam::phase::clone() const
{
    notImplemented("phase::clone() const");
    return autoPtr<phase>(NULL);
}

Foam::tmp<Foam::volScalarField> Foam::phase::mu
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return nuModel_->nu() * rho(p,T);
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
    nuModel_->correct();
}


bool Foam::phase::read(const dictionary& phaseDict)
{
    phaseDict_ = phaseDict;

    if (nuModel_->read(phaseDict_))
    {
        return true;
    }
    else
    {
        return false;
    }
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

// This calculates the sensible enthalpy of the phase. I do not think it is
// ever used
/*
Foam::tmp<volScalarField> Foam::phase::hs
(
    const volScalarField& T
) const
{
    Foam::Info<< "Calling phase::hs(T)" << Foam::endl;
    tmp<volScalarField> ths
    (
        new volScalarField
        (
            IOobject
            (
                "ths"+name_,
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("ths"+name_, dimEnergy/dimMass, 0.0)
        )
    );
    
    forAllConstIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {
        ths() += specieI().Y() * specieI().hs(T);
    }
    
    return ths;
}
*/

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
        tpsi() = W()/(dimensionedScalar("R", dimensionSet(1, 2, -2, -1, -1), 8314) * T);
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


// Construct the molar mass field. Should be zero if Yp == 0. If Yp == VSMALL
// though then W -> VSMALL
Foam::tmp<volScalarField> Foam::phase::W() const
{
    return Yp()/Np();
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



// Returns a mask of ones everywhere "in" the phase (to limit diffusion across
// the phase boundary)
Foam::tmp<Foam::volScalarField> Foam::phase::phaseMask() const
{    
    return pos(sharp(0.01) - 0.95);
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

/*
// Lennerd Jones diffusion coefficient calculation - from Swan
Foam::tmp<Foam::volScalarField> Foam::phase::Dk
(
    const word& specieKName
) const
{
    // Dij = sum((1-xk)/(xk/Dki) i ne k)
    tmp<volScalarField> trDk
    (
        new volScalarField
        (
            IOobject
            (
                "trDij"+name_,
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("rDij",dimTime/dimArea,SMALL)
        )
    );
    
    const subSpecie& specieK = *subSpecies_[specieKName];
    
    tmp<volScalarField> xk = specieK.Y() / (specieK.W() * Np())


    forAllConstIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {
        if( specieI().name() != specieKName )
        {
            trDk() += specieI().Y() / (specieI().W() * Np()) / Dij( specieI(), specieK );
        }
    }

    return (1 - xk()) / trDk;
}
*/


tmp<volVectorField> Foam::phase::calculateDs(scalar maskTol, bool allow)
{
    tmp<volVectorField> tUc
    (
        new volVectorField
        (
            IOobject
            (
                "tUc"+name_,
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedVector("Uc",dimVelocity,vector::zero)
        )
    );
        
    // DiM = sum((1-xk)/(xk/Dki) i ne k)  <-- TODO: VERIFY THIS!
    const phase& alpha = *this;
    tmp<volScalarField> tYp = Yp();
    
    forAllIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {
        tmp<volScalarField> tDen
        (
            new volScalarField
            (
                IOobject
                (
                    "tDen"+name_,
                    mesh().time().timeName(),
                    mesh()
                ),
                mesh(),
                dimensionedScalar("tDen",dimTime/dimArea,SMALL)
            )
        );
        
        if( allow )
        {
        
            tmp<volScalarField> xI = specieI().Y() / (specieI().W() * Np());
            
            forAllConstIter(PtrDictionary<subSpecie>, subSpecies_, specieJ)
            {
                if( specieJ().name() != specieI().name() )
                {
                    tDen() += specieI().Y() / (specieI().W() * Np()) / specieI().Dij( specieJ() );
                }
            }
            
            specieI().D() = (1.0 - xI) / tDen * pos(alpha - maskTol);
            
            const volScalarField& Yi = specieI().Y();
            tUc() += specieI().D()  * fvc::grad(Yi / (tYp()+SMALL));
        }
        else
        {
            specieI().D() = 0.0/tDen;
        }
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
void Foam::phase::solveSubSpecies
(
    const volScalarField& rhoTotal,
    const surfaceScalarField& rhoPhi,
    const volScalarField& p,
    const volScalarField& T,
    const phase& alphaLiquid,
    const volVectorField& uc,
    const tmp<fv::convectionScheme<scalar> >& mvConvection
)
{   
    word divScheme("div(rho*phi*alpha,Yi)");
    word divSchemeCorr("div(uc*Yi)");
    
    Info<< "Solving "<<subSpecies_.size()
        <<" subspecie(s) for phase: " 
        << name_ << endl;
        
    /*
        Approaches:
            next up - adding back turbulent mu (only in purely vapor cells, alphaLiquid < SMALL - see how much alphaL diffuses in current scheme)
            move DiM and uc into phase (it no longer needs to be in thermo)
            
            checkerboarding could be due to implicit evaporation source.
                this makes Sv and Sl inconsistent, and implicitly allows
                  condensation, but may be oscillatory. If current case still
                  checkboards, then try this avenue.
                  
                  Perhaps allowing condensation would be the best approach?
                  
              changing S_evap to explicit, checkerboarding does eventually go away tho
              
              Sp_Yi_evap uses Np, so it is changing as we loop through species. - Fix this
                fully explicit now so this isn't used, if it works, we can make it "smart implicit"
                adding condensation did not remove checkerboarding
                
            explicit, no Np issue still checkerboards!!
            
            trying it with no diffusion to see if that is the source of the checkerboarding
            
            currently running
                Dmask is based on alpha > 0.01
                    try a non-overlapping mask (alpha > 0.9) - boundedness much improved at first, maybe not later...
            
                let it run for awhile and look at
                    1. how it does
                    2. how a 0.99 mask would look (too diffuse?)
    */

    if (subSpecies_.size() > 1)
    {        
        tmp<volScalarField> Yp0 = Yp();
        
        forAllIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
        {
            volScalarField& Yi = specieI().Y();
            
            fvScalarMatrix YiEqn
            (
                fvm::ddt(rhoTotal, Yi)
              //+ fvm::div(rhoPhi, Yi, divScheme)
              + mvConvection->fvmDiv(rhoPhi, Yi)
              - fvc::laplacian(specieI().D()*rho(p,T), Yi/(Yp0()+SMALL)) //need to relax alpha field before this can work stably
              + fvc::div(uc*Yi*rho(p,T), divSchemeCorr)
             ==
                combustionPtr_->R(Yi)
              + Su_Yi_evap( alphaLiquid, specieI() )
              //- Sp_Yi_evap( alphaLiquid, specieI() )*Yi
              //- fvm::Sp( Sp_Yi_evap( alphaLiquid, specieI()), Yi)
            );

            YiEqn.relax();
            YiEqn.solve(mesh().solver("Yi"));
            
            Yi.max(0.0);
            Yi.min(1.0); //TODO: Try removing this
            
            Info<< "  Pre-coerce: " << Yi.name() << " min,max,avg = " 
                << Foam::min(Yi).value() <<  ", " << Foam::max(Yi).value() 
                << ", " << Yi.weightedAverage(mesh().V()).value() << endl;
        }
        
        
        // Coerce subspecies to host phase
        // This can only constrict (reduce) Ys, not expand (in Vapors)

/*
        tmp<volScalarField> Yp0 = Yp();
        tmp<volScalarField> YpMax = rho(p,T)/rhoTotal*sharp(0.0001);
        
        tmp<volScalarField> factor = YpMax / (Yp0+SMALL);
        
        if( name_ == "Vapor" )
        { //Do not expand vapors, but liquid phases must be expanded to fill host phase
            factor() = pos(factor() - 1.0) + neg(factor() - 1.0 + VSMALL)*factor();
        }
        
        forAllIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
        {
            volScalarField& Yi = specieI().Y();
            
            Yi *= factor();
            
            Yi.max(0.0);
            Yi.min(1.0);
            
            Info<< "  " << Yi.name() << " min,max,avg = " 
                << Foam::min(Yi).value() <<  ", " << Foam::max(Yi).value() 
                << ", " << Yi.weightedAverage(mesh().V()).value() << endl;
        }
        
        */
        
    }
    else
    {
        Info<<"  Setting Y for phase " << name() << endl;
        PtrDictionary<subSpecie>::iterator specieI = subSpecies_.begin();
        specieI().Y() = rho(p,T)/rhoTotal*sharp(0.0001);
    }
    
}



// ************************************************************************* //
