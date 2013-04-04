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
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Ypsum_"+name, dimless, 0.0)
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
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
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
    cellMask_.oldTime();
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
            //dimensionedScalar s("s",dimDensity,SMALL);
            //specieI().Yp() = specieI().Y()*cellMask_*(rhoAlpha_+otherRhoAlpha)/(rhoAlpha_+s);
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
        }
    }


/*
    forAllIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {   
        specieI().Yp() = specieI().Y() / (Yp() + SMALL);
        
        //Allow a little diffusion within masked region
        dimensionedScalar dA = Foam::pow(Foam::min(mesh().V()),2.0/3.0)/40.0;
    
        for(label i = 0; i < 5; i++)
        {
            specieI().Yp() += dA * fvc::laplacian(faceMask_, specieI().Yp());
        }
        specieI().Yp().min(1.0);
        specieI().Yp().max(0.0);
    }
    
    tmp<volScalarField> YpSum = Ypp();
    
    forAllIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    { 
        specieI().Yp() /= (YpSum() + SMALL);
        specieI().Yp().correctBoundaryConditions();
        specieI().Yp().oldTime();
    }*/
}

// Only applicable for liquid phase. Vapor phase will return 0 here.
// TODO: Have vapor specie do its sutherland calc here too
Foam::tmp<Foam::volScalarField> Foam::phase::mu
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    scalar muValue = (name_ == "Liquid") ? 1e-3 : 2e-5;
    
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
            dimensionedScalar("tmu", dimArea*dimDensity/dimTime, muValue)
        )
    );
    /*
    forAllConstIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {   
        if( specieI().hasNuModel() )
        {
            tmu() += specieI().nuModel().nu() * specieI().Yp() * rho(p,T);
        }
    }
    */
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
    rhoAlpha_ = *this * rho(p,T);// * cellMask_;
    
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
            
            tS_evap() += specieI().evapModel().rho_evap()*(Ru*T/(p*Wv) - 1/rho0);
        }
    }
    
    return tS_evap;
}

// Calculates the total mass source/sink due to evaporation.
Foam::tmp<Foam::volScalarField> Foam::phase::m_evap_sum() const
{
    tmp<volScalarField> tS_evap
    (
        new volScalarField
        (
            IOobject
            (
                "tm_evap",
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("m_evap", dimDensity/dimTime, 0.0)
        )
    );
    
    forAllConstIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {
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
                    if (specieLI().evapModel().vaporName() == specieI().Y().name())
                    {
                        tS_evap() += specieLI().evapModel().rho_evap();
                        break;
                    }
                }
            }
        }
        else
        {
            if (specieI().hasEvaporation())
            {
                tS_evap() -= specieI().evapModel().rho_evap();
            }
        }
    }
    
    return tS_evap;
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


Foam::Pair<Foam::tmp<Foam::volScalarField> > Foam::phase::pSuSp
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    Pair<tmp<volScalarField> > tSuSp
    (
        tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "tpSu",
                    mesh().time().timeName(),
                    mesh()
                ),
                mesh(),
                dimensionedScalar("pSu", dimless/dimTime, 0.0)
            )
        ),
        tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "tpSp",
                    mesh().time().timeName(),
                    mesh()
                ),
                mesh(),
                dimensionedScalar("pSp", dimless/dimTime/dimPressure, 0.0)
            )
        )
    );
    
    dimensionedScalar Ru("Ru",dimensionSet(1, 2, -2, -1, -1),8314);
    
    forAllConstIter(PtrDictionary<subSpecie>, subSpecies(), specieI)
    {
        dimensionedScalar rho0 = specieI().rho0();
        
        if (specieI().hasEvaporation() && rho0.value() > SMALL )
        {
            dimensionedScalar Wv = specieI().W();
            
            Pair<tmp<volScalarField> > tpSuSp = specieI().evapModel().pSuSp();
            
            
            //full source term linearization, including density relationship
            tmp<volScalarField> tdSdp = tpSuSp.second()()/rho0 - tpSuSp.first()()*Ru*T/(p*p*Wv);
            tmp<volScalarField> tSu = (tpSuSp.first()() - tpSuSp.second()()*p)*(Ru*T/(p*Wv) - 1/rho0) - tdSdp()*p;
            
            tSuSp.first()() += tSu;
            tSuSp.second()() -= tdSdp;
            
            
            //fully explicit treatment
            //tSuSp.first()() += (tpSuSp.first() - tpSuSp.second()*p)*(Ru*T/(p*Wv) - 1/rho0);
            //tSuSp.first()() += tpSuSp.first()*(Ru*T/(p*Wv) - 1/rho0);
            
        }
    }
        
    return tSuSp;
}


Foam::Pair<Foam::tmp<Foam::volScalarField> > Foam::phase::TSuSp() const
{
    Pair<tmp<volScalarField> > tSuSp
    (
        tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "tTSu",
                    mesh().time().timeName(),
                    mesh()
                ),
                mesh(),
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
                    mesh().time().timeName(),
                    mesh()
                ),
                mesh(),
                dimensionedScalar("TSp", dimPower/dimVolume/dimTemperature, 0.0)
            )
        )
    );
    
    forAllConstIter(PtrDictionary<subSpecie>, subSpecies(), specieI)
    {
        if( specieI().hasEvaporation() )
        {
            Pair<tmp<volScalarField> > tTSuSp = specieI().evapModel().TSuSp();
            tSuSp.first()() += tTSuSp.first();
            tSuSp.second()() += tTSuSp.second();
        }
    }
        
    return tSuSp;
}




Foam::Pair<Foam::tmp<Foam::volScalarField> > Foam::phase::YiSuSp
(
    const subSpecie& specieI
) const
{
    Pair<tmp<volScalarField> > tSuSp
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
                    tSuSp = specieLI().evapModel().YSuSp();
                    break;
                }
            }
        }
    }
    else
    {
        if (specieI.hasEvaporation())
        {
            //tSuSp = specieI.evapModel().YSuSp();
            tSuSp.first() = -specieI.evapModel().rho_evap();
        }
    }
    
    return tSuSp;
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
                    tSu_evap() += specieLI().evapModel().rho_evap();
                    break;
                }
            }
        }
    }
    else
    {
        if (specieI.hasEvaporation())
        {
            tSu_evap() -= specieI.evapModel().rho_evap();
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
            dimensionedScalar("sigma", dimensionSet(1, 0, -2, 0, 0), 0.07)
        )
    );
    
    /*
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

    return tsigma/tn;*/
    return tsigma;
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


Foam::tmp<volVectorField> Foam::phase::URecoil() const
{
    tmp<volVectorField> tUr
    (
        new volVectorField
        (
            IOobject
            (
                "tUr"+name_,
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedVector("tUr"+name_, dimVelocity, vector::zero)
        )
    );

    forAllConstIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
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
                if (specieLI().evapModel().vaporName() == specieI().Y().name())
                {
                    tUr() += specieLI().evapModel().U_evap();
                    break;
                }
            }
        }
    }
    
    Info<<"Max recoil velocity = " << Foam::max(Foam::mag(tUr())).value() << endl;

    return tUr;
}

Foam::tmp<volScalarField> Foam::phase::pRecoil() const
{
    tmp<volScalarField> tpR
    (
        new volScalarField
        (
            IOobject
            (
                "tpR"+name_,
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("tpR"+name_, dimPressure, 0.0)
        )
    );

    forAllConstIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
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
                if (specieLI().evapModel().vaporName() == specieI().Y().name())
                {
                    tpR() += specieLI().evapModel().pRecoil();
                    break;
                }
            }
        }
    }
    
    return tpR;
}

// Volumetric sink term from evaporation = rho_evap/rhoL (only applicable to
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
            tSv() += specieI().evapModel().rho_evap() / specieI().rho0();
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
        
        /*tmp<volScalarField> Yp_ = Ypp();
        tmp<volScalarField> Yvoid = 0.0001*neg(Yp_()-0.05);
        tmp<volScalarField> den = Yvoid()/rhoBase;
        
        forAllConstIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
        {
            den() += specieI().Yp() / specieI().rho0();
        }
        
        return (Yvoid + Yp_)/den;*/
        trho() = rhoBase;
        return trho;
        
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
        dimensionedScalar W0("W0",dimMass/dimMoles, 28);
        tpsi() = W0 / (Ru * T);
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
    
    surfaceScalarField DbyDelta // m/s diffusion velocity
    (
        mesh().surfaceInterpolation::deltaCoeffs()
      * fvc::interpolate(DbyRho)
    );
    
    scalar DiNum = gMax(DbyDelta.internalField())*mesh().time().deltaTValue();
    
    Info<< "Max DiNum = " << DiNum << endl;
        
    Info<< "Solving "<<subSpecies_.size()
        <<" subspecie(s) for phase: " 
        << name_ << endl;
        
    //solve ddt(rhoAlpha) + div(rhoPhiAlpha) == 0
    /*scalarField& rhoAlphaIf = rhoAlpha_;
    const scalarField& rhoAlpha0 = rhoAlpha_.oldTime();
    const scalar deltaT = mesh().time().deltaTValue();

    rhoAlphaIf = 0.0;
    surfaceScalarField rhoPhiAlphaMasked = rhoPhiAlpha_*faceMask_;
    fvc::surfaceIntegrate(rhoAlphaIf, rhoPhiAlphaMasked);
    
    volScalarField Su = m_evap_sum() * cellMask_;
    
    forAll(rhoAlphaIf, cellI)
    {
        if( cellMask_[cellI] < SMALL && mag(Su[cellI]) > SMALL )
        {
            Info<<"WARNING: evaporation outside mask region: "<<Su[cellI]
                <<" at "<<mesh().C()[cellI]<<endl;
        }
        
        //new amount = old amount + (source - outflux) * dt
        rhoAlphaIf[cellI] = rhoAlpha0[cellI] + (Su[cellI] - rhoAlphaIf[cellI])*deltaT;
        

    }
    rhoAlpha_.correctBoundaryConditions();*/
    
    //Arbitrary diagonal term for cells outside normal region
    volScalarField Sp = (1.0 - cellMask_)*rho(p,T)/mesh().time().deltaT();
    
    //Zero mass outside region
    //rhoAlpha_ *= cellMask_;
    rhoPhiAlpha_ *= faceMask_;
    
   
    // Loop through phase's subspecies
    //scalar dTsrc = mesh().time().deltaTValue();
    
    forAllIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {
        Info<<"Solving specie " << specieI().Y().name() << endl;
        
        volScalarField& Yi = specieI().Yp();
        
        //tmp<volScalarField> SuEvap = Su_Yi_evap( specieI() );
        
        Pair<tmp<volScalarField> > YSuSp = YiSuSp( specieI() );
        
        //dimensionedScalar ss("ss",dimDensity/dimTime,SMALL);
        //scalar mindT = Foam::min(0.02 * rho(p,T) / (Foam::mag(SuEvap())+ss)).value();
        //dTsrc = (mindT < dTsrc) ? mindT : dTsrc;
        
        const rhoChemistryModel& chemistry = combustionPtr_->pChemistry();
        const volScalarField& kappa = mesh().lookupObject<volScalarField>("PaSR::kappa");
        
        tmp<volScalarField> R = kappa * chemistry.RR( specieI().idx() );
        
        fvScalarMatrix YiEqn
        (
            fvm::ddt(rhoAlpha_, Yi)
          + fvm::div(rhoPhiAlpha_, Yi, divScheme)
          - fvm::Sp(fvc::ddt(rhoAlpha_) + fvc::div(rhoPhiAlpha_) - m_evap_sum(), Yi)
          - fvm::laplacian(D_, Yi)
         ==
            R()
          + YSuSp.first()
          - fvm::Sp(YSuSp.second(), Yi)
          - fvm::Sp(Sp, Yi)
        );
        
        Info<<"YiEqn.A() min = " << Foam::min(YiEqn.A()) << endl;
        
        if( Foam::min(YiEqn.A()).value() < SMALL )
        {
            tmp<volScalarField> A = YiEqn.A();
            tmp<volScalarField> ddtdiv = fvc::ddt(rhoAlpha_) + fvc::div(rhoPhiAlpha_) - m_evap_sum();
            tmp<volScalarField> div = fvc::div(rhoPhiAlpha_);
            tmp<volScalarField> ddt = fvc::ddt(rhoAlpha_);
            tmp<volScalarField> avgFaceMask = fvc::average(faceMask_);
            tmp<volScalarField> avgD = fvc::average(D_);
            
            forAll(A(), cellI)
            {
                if( A()[cellI] < SMALL )
                {
                    
                    Info<<"A = " << A()[cellI] << " where" << endl;
                    Info<<"  rhoAlpha = " << rhoAlpha_[cellI] << endl;
                    Info<<"  rhoAlphaOld = " << rhoAlpha_.oldTime()[cellI] << endl;
                    Info<<"  Sp = " << Sp[cellI] << endl;
                    Info<<"  cellMask = " << cellMask_[cellI] << endl;
                    Info<<"  cellMaskOld = " << cellMask_.oldTime()[cellI] << endl;
                    Info<<"  alpha = " << this->operator[](cellI) << endl;
                    Info<<"  alphaOld = " << this->oldTime()[cellI] << endl;
                    Info<<"  ddtdiv = " << ddtdiv()[cellI] << endl;
                    Info<<"  ddt = " << ddt()[cellI] << endl;
                    Info<<"  div = " << div()[cellI] << endl;
                    Info<<"  avgFM = " << avgFaceMask()[cellI] << endl;
                    Info<<"  avgD = " << avgD()[cellI] << endl;
                    Info<<"  Y = " << Yi[cellI] << endl;
                }
            
            }
        
        }

        //Info<<"Doing solve, min diag = " << Foam::min(YiEqn.A()) << endl;
        YiEqn.relax();
        YiEqn.solve(mesh().solver("Yi"));
        
        Yi.max(0.0);
        Yi.min(1.0);
        
        Info<< Yi.name() << " min,max,avg = " 
            << Foam::min(Yi).value() <<  ", " << Foam::max(Yi).value() 
            << ", " << Yi.weightedAverage(mesh().V()).value() << endl;  
    }
    
    
    /*tmp<volScalarField> Ysum = Ypp() + SMALL;
    
    forAllIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {
        volScalarField& Yi = specieI().Yp();
        
        Yi *= cellMask_/Ysum();
        
        Info<< Yi.name() << " min,max,avg = " 
            << Foam::min(Yi).value() <<  ", " << Foam::max(Yi).value() 
            << ", " << Yi.weightedAverage(mesh().V()).value() << endl;  
    }*/
    
    //Info<< "Min source dT = " << dTsrc << endl;
     
    return DiNum;
}


// normal vector pointing outward from this phase
Foam::tmp<Foam::volVectorField> Foam::phase::n() const
{
    dimensionedScalar deltaN
    (
        "deltaN",
        1e-8/Foam::pow(Foam::average(mesh().V()), 1.0/3.0)
    );
    
    tmp<volVectorField> nOut = -fvc::grad(*this);
    nOut() /= (mag(nOut()) + deltaN);
    
    return nOut;
}

// Update the phase masks based on alpha and the evaporation area
void Foam::phase::setPhaseMasks
(
    scalar maskTol,
    const volScalarField& p,
    const volScalarField& T
)
{
    volScalarField& alpha = *this;
    
    //TODO: Make this independent of specie type "N2O4"
    //const volScalarField& area = mesh().lookupObject<volScalarField>("area_N2O4");
       
       
    const volScalarField* areaPtr = NULL;
    if( name_ == "Vapor" )
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
                areaPtr = &(specieLI().evapModel().area());
                break;
            }
        }
    }
    else
    {
        forAllConstIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
        {
            if (specieI().hasEvaporation())
            {
                areaPtr = &(specieI().evapModel().area());
                break;
            }
        }
    }
    
    if( areaPtr == NULL )
    {
        Info<< "WARNING: NULL area pointer"<<endl;
    }
    
    const volScalarField& area = *areaPtr;
       
       
       
       
    dimensionedScalar one("one",dimLength,1.0);
    
    //tmp<volScalarField> oldMask = cellMask_;
    cellMask_ = pos(alpha - maskTol + area*one - SMALL);
    faceMask_ = pos(fvc::interpolate(cellMask_) - 0.95);
    
    //apply the current-time mask to the old-time rhoAlpha
    // ddt(rhoAlpha) = ddt(rho*alpha*mask) 
    //               = ddt(mask)*rho*alpha + mask*ddt(rho*alpha)
    //
    // The ddt(mask) term is either 0 or 1/dt, and this extremely large change
    // in magnitude makes the mask transition unstable. The time rate of change
    // of the mask is irrelevant to the computation, we only really want the
    // mask*ddt(rho*alpha) term, so we remove this by retro-actively applying
    // the mask to oldTime.
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
    //rhoAlpha_.oldTime() = cellMask_ * alpha.oldTime() * rho(p.oldTime(), T.oldTime());
    
    //Perform a local averaging to "fill in" newly un-masked cells
    //forAllIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    //{ 
        //tmp<volScalarField> tYpf = fvc::average(fvc::interpolate(specieI().Yp()));
        
    /*if( name_ == "Liquid")
    {
        forAll(cellMask_, cellI)
        {
            if( cellMask_[cellI] > 0.5 && cellMask_.oldTime()[cellI] < 0.5 )
            {
                //specieI().Yp()[cellI] += tYpf()[cellI];
                rhoAlpha_.oldTime()[cellI] = rhoAlpha_[cellI]; 
            }
        }
    }*/
    //}
}

void Foam::phase::updateGlobalYs
(
    const volScalarField& rhoAlphaOther
)
{
    forAllIter(PtrDictionary<subSpecie>, subSpecies_, specieI)
    {   
        specieI().Y() = specieI().Yp()*rhoAlpha_/(rhoAlpha_ + rhoAlphaOther);
        specieI().Y().max(0.0);
        specieI().Y().correctBoundaryConditions();
    }
    
    Ypsum_ = Ypp();
}


// ************************************************************************* //
