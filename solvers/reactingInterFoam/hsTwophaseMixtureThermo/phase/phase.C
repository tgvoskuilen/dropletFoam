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
#include "mixturePhaseChangeModel.H"

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
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rhoAlpha_"+name, dimDensity, 1.0)
    ),
    Ypsum_
    (
        IOobject
        (
            "Ypsum_"+name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("Ypsum_"+name, dimless, 0.0)
    ),
    otherPhase_(NULL),
    combustionPtr_(NULL),
    species_(species),
    speciesData_(speciesData),
    Sc_
    (
        phaseDict_.lookupOrDefault
        (
            "Sc",
            dimensionedScalar("Sc",dimless,1.0)
        )
    ),
    subSpecies_
    (
        phaseDict_.lookup("subspecies"),
        subSpecie::iNew(mesh, species, speciesData, Sc_)
    ),
    cellMask_
    (
        IOobject
        (
            "cellMask_"+name,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("cellMask_"+name, dimless, 1.0)
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
        dimensionedScalar("faceMask_"+name, dimless, 1.0)
    )
{  
    this->oldTime();
    cellMask_.oldTime();
    
    Info<< "Created phase " << name << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::phase::setOtherPhase
(
    const phase* other
)
{
    otherPhase_ = other;
}
        
        
Foam::autoPtr<Foam::phase> Foam::phase::clone() const
{
    notImplemented("phase::clone() const");
    return autoPtr<phase>(NULL);
}

void Foam::phase::setSpecies( const volScalarField& otherRhoAlpha )
{    
    bool isIC = false;
    
    forAllIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {
        //Y == Yp and Yp != 0
        if( Foam::max(Foam::mag(specieI().Yp() - specieI().Y())).value() < 1e-4 
            && Foam::max(Foam::mag(specieI().Yp())).value() > 0.01 )
        {
            Info<<"Initializing Yp_" << specieI().name() << " from Y" << endl;
            isIC = true;
            
            //only set Yp if Y == Yp (only the case if Yp not read from file)
            specieI().Yp() = specieI().Y() / (Yp() + SMALL);
            
            //Allow a little diffusion within masked region
            dimensionedScalar dA = Foam::pow(Foam::min(mesh().V()),2.0/3.0)/40.0;
        
            for(label i = 0; i < 5; i++)
            {
                specieI().Yp() += dA*fvc::laplacian(faceMask_, specieI().Yp());
            }
            specieI().Yp().min(1.0);
            specieI().Yp().max(0.0);
        }
    }
    
    if( isIC )
    {
        tmp<volScalarField> YpSum = Ypp();
        
        forAllIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
        { 
            specieI().Yp() *= cellMask_/(YpSum() + SMALL);
            specieI().Yp().correctBoundaryConditions();
            specieI().Yp().oldTime();
            
            Info<< "Min,Max values of " << specieI().name() << " = " 
                << Foam::min(specieI().Yp()).value() << ", " 
                << Foam::max(specieI().Yp()).value() << endl;
        }
    }


    //Initialize the phase density field
    const volScalarField& p = mesh().lookupObject<volScalarField>("p");
    const volScalarField& T = mesh().lookupObject<volScalarField>("T");
    
    rhoAlpha_ = sharp(0.0)*rho(p,T)*cellMask_;
    
    rhoAlpha_.oldTime();
}

// Calculate viscosity
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
            // Use viscosity transportModel (liquids)
            tmu() += specieI().nuModel().nu() * specieI().Yp() * rho(p,T);
        }
        else
        {
            // Use Sutherland model (vapors)
            tmu() += specieI().mu(T) * specieI().Yp();
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
    rhoAlpha_ = sharp(0.0) * rho(p,T) * cellMask_;
    rhoAlpha_.correctBoundaryConditions();
        
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



// Calculate the mole fraction of a named specie
Foam::tmp<Foam::volScalarField> Foam::phase::x(const word& specie) const
{
    tmp<volScalarField> tx
    (
        new volScalarField
        (
            IOobject
            (
                "tx",
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("x", dimless, 0.0)
        )
    );

    forAllConstIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {
        if( specieI().name() == specie )
        {
            const volScalarField& Yi = specieI().Yp();
            dimensionedScalar Wi = specieI().W();
            tx() = Yi / (Wi * Npp());
            break;
        }
    }

    return tx;
}

// Calculate x/Y of a named specie
Foam::tmp<Foam::volScalarField> Foam::phase::xByY(const word& specie) const
{
    tmp<volScalarField> txByY
    (
        new volScalarField
        (
            IOobject
            (
                "txByY",
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("xByY", dimless, 0.0)
        )
    );

    forAllConstIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {
        if( specieI().name() == specie )
        {
            dimensionedScalar Wi = specieI().W();
            txByY() = 1.0 / (Wi * Npp());
            break;
        }
    }

    return txByY;
}



Foam::Pair<Foam::tmp<Foam::volScalarField> > Foam::phase::YiSuSp
(
    const subSpecie& specieI,
    const PtrDictionary<mixturePhaseChangeModel>& phaseChangeModels
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
                    mesh().time().timeName(),
                    mesh()
                ),
                mesh(),
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
                    mesh().time().timeName(),
                    mesh()
                ),
                mesh(),
                dimensionedScalar("YSp", dimDensity/dimTime, 0.0)
            )
        )
    );
    
    forAllConstIter
    (
        PtrDictionary<mixturePhaseChangeModel>, 
        phaseChangeModels, 
        pcmI
    )
    {
        if( pcmI().hasSpecie( specieI.name() ) )
        {
            Pair<tmp<volScalarField> > pcmYiSuSp = pcmI().YSuSp(specieI.name());
            tYSuSp.first()() += pcmYiSuSp.first();
            tYSuSp.second()() += pcmYiSuSp.second();
        }
    }
    
    return tYSuSp;
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
        
        tmp<volScalarField> Yp_ = Ypp();
        tmp<volScalarField> Yvoid = 0.0001*neg(Yp_()-0.05);
        tmp<volScalarField> den = Yvoid()/rhoBase;
        
        forAllConstIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
        {
            den() += specieI().Yp() / specieI().rho0();
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

Foam::tmp<volScalarField> Foam::phase::Npp() const
{
    //Npp = sum(Ypi/Wi)
    tmp<volScalarField> tNpp
    (
        new volScalarField
        (
            IOobject
            (
                "tNpp"+name_,
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("tNpp"+name_, dimensionSet(-1,0,0,0,1), SMALL)
        )
    );
    
    forAllConstIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {
        tNpp() += specieI().Yp() / specieI().W();
    }
    
    return tNpp;
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
    tmp<volScalarField> Wother = otherPhase_->Ypp() / otherPhase_->Npp();
    
    tmp<volScalarField> Yp_ = Ypp();
    tmp<volScalarField> Yvoid = 0.0001*neg(Yp_()-0.05);
    tmp<volScalarField> den = Yvoid()/(Wother + ws);
    
    forAllConstIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {
        den() += specieI().Yp() / specieI().W();
    }
    
    return (Yvoid + Yp_)/den;
}



// Calculate the phase thermal conductivity using the method from Harvazinski
Foam::tmp<Foam::volScalarField> Foam::phase::kappa
(
    const volScalarField& T
) const
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
        if( name_ == "Vapor" )
        {
            tk1() += specieI().Y() / (specieI().W() * Np()) * specieI().kappa(T);
            tk2() += specieI().Y() / (specieI().W() * Np()) / specieI().kappa(T);
        }
        else
        {
            tk1() += specieI().Y() / (specieI().W() * Np()) * specieI().kappaL();
            tk2() += specieI().Y() / (specieI().W() * Np()) / specieI().kappaL();
        }
    }
    
    tk1() += pos(Yp() - 1e-3) / tk2;

    return 0.5*tk1;
}




void Foam::phase::calculateDs
(
    const volScalarField& mut,
    const volScalarField& p,
    const volScalarField& T
)
{
    tmp<volScalarField> trho = rho(p,T);
    
    forAllIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {
        specieI().calculateDs(mut, trho(), T);
    }
}



// Calculate phase specific heat using subspecie Cp model, which calls the
// underlying thermo model for each specie (Janaf table)
// Where the phase is trace (Yp < 1e-4) the Cv is just the average of the
// phase's subspecies, otherwise it is weighted by mass fraction
Foam::tmp<Foam::volScalarField> Foam::phase::Cp
(
    const volScalarField& T
) const
{
    tmp<volScalarField> tCp
    (
        new volScalarField
        (
            IOobject
            (
                "tCp"+name_,
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("Cp",dimEnergy/dimMass/dimTemperature,SMALL)
        )
    );
    
    label ns = subSpecies_.size();
    
    tmp<volScalarField> tYp = Yp();
    tmp<volScalarField> mask1 = pos(tYp() - 1e-4);
    tmp<volScalarField> mask2 = (1 - mask1()) / scalar(ns);
    mask1() /= (tYp + SMALL);
    
    forAllConstIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {
        tCp() += ( mask1()*specieI().Y() + mask2() )*specieI().Cp(T);
    }
    
    tCp().correctBoundaryConditions();
    
    return tCp;
}

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
    scalar tolL,
    scalar tolH //defaults to -1, so same as tolL
) const
{
    if( tolH < 0.0 )
    {
        tolH = tolL;
    }
    
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
    scalar Cpc = tolL+tolH;
    s = (Foam::min(Foam::max(s, tolL),1.0-tolH) - tolL)/(1.0-Cpc);
    
    s.correctBoundaryConditions();

    return ts;
}

        
// Solve the mass fraction governing equation (YEqn) for each subspecie of
// this phase.
Foam::tmp<Foam::surfaceScalarField> Foam::phase::solveSubSpecies
(
    const volScalarField& p,
    const volScalarField& T,
    const PtrDictionary<mixturePhaseChangeModel>& phaseChangeModels
)
{   
    word divScheme("div(rho*phi*alpha,Yi)");
            
    Info<< "Solving "<<subSpecies_.size()
        <<" subspecie(s) for phase: " 
        << name_ << endl;
        
    // Add up the total mass generation for this phase due to phase change
    volScalarField mdot_phase
    (
        IOobject
        (
            "m_pc",
            mesh().time().timeName(),
            mesh()
        ),
        mesh(),
        dimensionedScalar("m_pc",dimDensity/dimTime,0.0)
    );

    forAllConstIter
    (
        PtrDictionary<mixturePhaseChangeModel>, 
        phaseChangeModels, 
        pcmI
    )
    {
        mdot_phase += pcmI().mdot(name_);
    }
        
           
    //Arbitrary diagonal term for cells outside normal region
    volScalarField Sp = (1.0 - cellMask_)*rho(p,T)/mesh().time().deltaT();
    
    // Create a container to add up diffusion-driven energy flux
    tmp<surfaceScalarField> tDgradYCp
    (
        new surfaceScalarField
        (
            IOobject
            (
                "DgradY",
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("DgradY",dimPower/dimTemperature,0.0)
        )
    );
    surfaceScalarField& DgradYCp = tDgradYCp();
    
    // Loop over all subspecies in this phase
    forAllIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {
        Info<<"Solving specie " << specieI().Y().name() << endl;
        
        volScalarField& Yi = specieI().Yp();
        const surfaceScalarField& Di = specieI().D();
        
        // Generate source term pair from phase change
        Pair<tmp<volScalarField> > YSuSp = YiSuSp( specieI(), phaseChangeModels );
        
        // Generate reaction-based source term
        const rhoChemistryModel& chemistry = combustionPtr_->pChemistry();
        const volScalarField& kappa = mesh().lookupObject<volScalarField>("PaSR::kappa");
        tmp<volScalarField> R = kappa * chemistry.RR( specieI().idx() );
        
        // Build specie governing equation
        fvScalarMatrix YiEqn
        (
            fvm::ddt(rhoAlpha_, Yi)
          + fvm::div(rhoPhiAlpha_, Yi, divScheme)
          - fvm::Sp(fvc::ddt(rhoAlpha_) + fvc::div(rhoPhiAlpha_) - mdot_phase, Yi)
          - fvm::laplacian(Di*faceMask_, Yi)
         ==
            R()
          + YSuSp.first()
          - fvm::SuSp(YSuSp.second(), Yi)
          - fvm::Sp(Sp, Yi)
        );
                
        // Solve specie mass fraction equation
        YiEqn.relax();
        YiEqn.solve(mesh().solver("Yi"));
        
        // Add up diffusion-driven energy flux
        DgradYCp += Di * faceMask_ * fvc::snGrad(Yi) * mesh().magSf() * fvc::interpolate(specieI().Cv(T));
        
        // A possibly more consistent way?
        //Efluxp += YiEqn.flux() * fvc::interpolate(specieI().Cv(T));
        
        Yi.max(0.0);
        Yi.min(1.0);
        
        Info<< Yi.name() << " min,max,avg = " 
            << Foam::min(Yi).value() <<  ", " << Foam::max(Yi).value() 
            << ", " << Yi.weightedAverage(mesh().V()).value() << endl;  
    }
     
    return tDgradYCp;
}



// Update the phase masks based on alpha and the evaporation area
void Foam::phase::setPhaseMasks
(
    scalar maskTol,
    const volScalarField& p,
    const volScalarField& T,
    const PtrDictionary<mixturePhaseChangeModel>& phaseChangeModels,
    const volScalarField& area
)
{    
    volScalarField& alpha = *this;
    
    // Add up the total mass source term for this phase due to phase change
    volScalarField mdot_phase
    (
        IOobject
        (
            "m_pc",
            mesh().time().timeName(),
            mesh()
        ),
        mesh(),
        dimensionedScalar("m_pc",dimDensity/dimTime,0.0)
    );
        
    forAllConstIter
    (
        PtrDictionary<mixturePhaseChangeModel>, 
        phaseChangeModels, 
        pcmI
    )
    {
        mdot_phase += pcmI().mdot(name_);
    }
    
    dimensionedScalar oneM("oneM",dimTime/dimDensity,1.0);
    dimensionedScalar oneA("oneA",dimLength,1.0);
    
    // Problem Corner Cases:
    //  1. alpha = 0, mdot = 0, but area > 0 (solution to Y equation will have 0 diagonal)
    //  2. alpha = 1e-8, alphaOld = 0, div = 0, laplacian = 0, area > 0
    //
    //  area condition seems too problematic
    cellMask_ = pos(alpha - maskTol + mag(mdot_phase)*oneM - SMALL);
    faceMask_ = pos(fvc::interpolate(cellMask_) - 0.95);
    
    //apply the current-time mask to the old-time rhoAlpha
    // ddt(rhoAlpha) = ddt(rho*alpha*mask) 
    //               = ddt(mask)*rho*alpha + mask*ddt(rho*alpha)
    //
    // The ddt(mask) term is either 0 or 1/dt, and this extremely large change
    // in magnitude makes the mask transition unstable. The time rate of change
    // of the mask is irrelevant to the computation, we only really want the
    // mask*ddt(rho*alpha) term, so we remove this by retro-actively applying
    // the mask to oldTime so ddt(mask) is always 0
    //
    tmp<volScalarField> ddtM = fvc::ddt(cellMask_);
    tmp<volScalarField> rhoOld = rho(p.oldTime(), T.oldTime());
    forAll(ddtM(), cellI)
    {
        if( ddtM()[cellI] > 10. )
        {
            rhoAlpha_.oldTime()[cellI] = cellMask_[cellI]
                 * alpha.oldTime()[cellI] * rhoOld()[cellI];
        }
    }
    
    cellMask_.correctBoundaryConditions();
    
    rhoPhiAlpha_ *= faceMask_;
    
}

void Foam::phase::updateGlobalYs
(
    const volScalarField& myRhoAlpha,
    const volScalarField& otherRhoAlpha
)
{
    dimensionedScalar s("small",dimDensity,SMALL);
    
    forAllIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {   
        specieI().Y() = specieI().Yp()*myRhoAlpha/(myRhoAlpha + otherRhoAlpha + s);
        specieI().Y().max(0.0);
        specieI().Y().correctBoundaryConditions();
    }
    
    Ypsum_ = Ypp();
}


// ************************************************************************* //
