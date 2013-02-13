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

#include "hsTwophaseMixtureThermo.H"
#include "fvMesh.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class MixtureType>
void Foam::hsTwophaseMixtureThermo<MixtureType>::calculate()
{
    const scalarField& hsCells = hs_.internalField();

    scalarField& TCells = T_.internalField();
    //scalarField& alphaCells = alpha_.internalField();
    
    Info << "Calculating mixture temperature from hs" << endl;

    forAll(TCells, celli)
    {
        // references multiComponentMixture.C
        const typename MixtureType::thermoType& mixture =
            this->cellMixture(celli);

        TCells[celli] = mixture.THs(hsCells[celli], TCells[celli]);
        //alphaCells[celli] = mixture.alpha(TCells[celli]);
    }

    forAll(T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pT = T_.boundaryField()[patchi];
        fvPatchScalarField& phs = hs_.boundaryField()[patchi];
        //fvPatchScalarField& palpha = alpha_.boundaryField()[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture =
                    this->patchFaceMixture(patchi, facei);

                phs[facei] = mixture.Hs(pT[facei]);

                //palpha[facei] = mixture.alpha(pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture =
                    this->patchFaceMixture(patchi, facei);

                pT[facei] = mixture.THs(phs[facei], pT[facei]);

                //palpha[facei] = mixture.alpha(pT[facei]);
            }
        }
    }
    
    
    calculateMu();
    calculateRho();
    calculatePsi();
}




// The internally stored "mu_" is used by turbulence and should be correct
template<class MixtureType>
void Foam::hsTwophaseMixtureThermo<MixtureType>::calculateMu()
{    
    vapor_.correct();
    liquid_.correct();
    
    mu_ = vapor_ * vapor_.mu() + liquid_ * liquid_.mu();
    mu_.correctBoundaryConditions();
}



template<class MixtureType>
void Foam::hsTwophaseMixtureThermo<MixtureType>::calculatePsi()
{    
    psi_ = vapor_*vapor_.psi();
    psi_.correctBoundaryConditions();
}


template<class MixtureType>
void Foam::hsTwophaseMixtureThermo<MixtureType>::calcEvaporation()
{
    Foam::Info << "Calculating evaporation rates" << Foam::endl;

    tmp<volScalarField> tAbyV
    (
        new volScalarField
        (
            IOobject
            (
                "tAbyV",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("tAbyV", dimless/dimLength, 0.0)
        )
    );
    
    tAbyV().internalField() = interface_.area().internalField() / mesh_.V();
    tAbyV().correctBoundaryConditions();

    forAllIter(PtrDictionary<subSpecie>, liquid_.subSpecies(), specieI)
    {
        if (specieI().hasEvaporation())
        {
            specieI().evapModel().calculate( tAbyV() );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MixtureType>
Foam::hsTwophaseMixtureThermo<MixtureType>::hsTwophaseMixtureThermo
(
    const fvMesh& mesh
)
:
    hsReactionThermo(mesh),
    MixtureType(*this, mesh),
    mesh_(mesh),
    combustionPtr_(NULL),
    vapor_
    (
        "Vapor",
        subDict("Vapor"),
        mesh,
        this->Y(),
        this->speciesData()
    ),
    liquid_
    (
        "Liquid",
        subDict("Liquid"),
        mesh,
        this->Y(),
        this->speciesData()
    ),
    phi_( mesh_.lookupObject<surfaceScalarField>("phi") ),
    U_( mesh_.lookupObject<volVectorField>("U") ),
    rhoPhi_
    (
        IOobject
        (
            "rhoPhi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rhoPhi", dimMass/dimTime, 0.0)
    ),
    sigma_
    (
        IOobject
        (
            "sigma",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("sigma", dimensionSet(1, 0, -2, 0, 0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    YSum_
    (
        IOobject
        (
            "Ysum",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Ysum", dimless, 0.0)
    ),
    deltaN_
    (
        "deltaN",
        1e-8/pow(average(mesh.V()), 1.0/3.0)
    ),
    phaseRelaxTime_(readScalar(lookup("relaxTime"))),
    interface_
    (
        liquid_,
        vapor_,
        U_
    )
{

    // Create convection field table
    /*forAll(this->Y(), i)
    {
        fields_.add(this->Y()[i]);
    }*/

    //Set evaporation models
    Info<< "Setting evaporation models" << endl;

    forAllIter(PtrDictionary<subSpecie>, liquid_.subSpecies(), specieI)
    {
        if(specieI().hasEvaporation())
        {
            specieI().evapPtr() = evaporationModel::New
            (
                specieI().dict(),
                p_,
                T_,
                liquid_,
                vapor_
            );
        }
    }
    
    calculateRho();
    rho_.oldTime();   
        
    setHs(T_);

    calculate();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class MixtureType>
Foam::hsTwophaseMixtureThermo<MixtureType>::~hsTwophaseMixtureThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Update the stored density field (rho_)
template<class MixtureType>
void Foam::hsTwophaseMixtureThermo<MixtureType>::calculateRho()
{
    liquid_.correct();
    vapor_.correct();
    
    rho_ = liquid_*liquid_.rho() + vapor_*vapor_.rho();
    rho_.correctBoundaryConditions();
    
    if( Foam::min(rho_).value() < SMALL )
    {
        FatalErrorIn
        (
            "multiphaseReactingMixture::calculateRho()"
        )   << "Encountered a zero mixture density"
            << "\n    Check that phase boundary conditions have been applied"
            << "\n    and that a nonzero alpha is specified on each boundary."
            << exit(FatalError);
    }
}


template<class MixtureType> 
Foam::tmp<Foam::volScalarField> 
Foam::hsTwophaseMixtureThermo<MixtureType>::W() const
{
    tmp<volScalarField> tden
    (
        new volScalarField
        (
            IOobject
            (
                "tden",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("tden", dimMoles/dimMass, 0.0)
        )
    );

    forAllConstIter(PtrDictionary<subSpecie>, vapor_.subSpecies(), specieI)
    {
        tden() += specieI().Y() / specieI().W();
    }
    
    forAllConstIter(PtrDictionary<subSpecie>, liquid_.subSpecies(), specieI)
    {
        tden() += specieI().Y() / specieI().W();
    }

    return 1.0/tden();
}


// Cell centered mixture viscosity
template<class MixtureType> 
Foam::tmp<Foam::volScalarField> 
Foam::hsTwophaseMixtureThermo<MixtureType>::muV() const
{
    tmp<volScalarField> tmu
    (
        new volScalarField
        (
            IOobject
            (
                "tmu",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("tmu", dimMass/dimTime/dimLength, 0.0)
        )
    );


    scalarField& muCells = tmu().internalField();
    const scalarField& TCells = T_.internalField();

    forAll(TCells, celli)
    {
        muCells[celli] = this->cellMixture(celli).mu( TCells[celli] );
    }

    forAll(T_.boundaryField(), patchi)
    {
        tmu().boundaryField()[patchi] = mu(T_.boundaryField()[patchi], patchi);
    }

    return tmu;
}

template<class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::hsTwophaseMixtureThermo<MixtureType>::mu
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tmu(new scalarField(T.size()));

    scalarField& mu = tmu();

    forAll(T, facei)
    {
        mu[facei] = this->patchFaceMixture(patchi, facei).mu(T[facei]);
    }

    return tmu;
}



// Face-centered viscosity
template<class MixtureType> 
Foam::tmp<Foam::surfaceScalarField> 
Foam::hsTwophaseMixtureThermo<MixtureType>::muf() const
{
    tmp<surfaceScalarField> tmuf =
        fvc::interpolate(vapor_)*fvc::interpolate(muV())
      + fvc::interpolate(liquid_)*fvc::interpolate(liquid_.mu(p_,T_));

    return tmuf;
}


// Total enthalpy source/sink, including reactions and evaporation
template<class MixtureType> 
Foam::tmp<Foam::volScalarField> 
Foam::hsTwophaseMixtureThermo<MixtureType>::Sh() const
{
    return combustionPtr_->Sh() + liquid_.Sh_evap();
}



// Used for alphaCourant number calculations
template<class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::hsTwophaseMixtureThermo<MixtureType>::nearInterface
(
    double lower, 
    double upper
) const
{
    tmp<volScalarField> tnearInt
    (
        new volScalarField
        (
            IOobject
            (
                "nearInterface",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("nearInterface", dimless, 0.0)
        )
    );


    tnearInt() = max(tnearInt(), pos(liquid_-lower)*pos(upper-liquid_));
    tnearInt() = max(tnearInt(), pos(vapor_-lower)*pos(upper-vapor_));

    return tnearInt;
}



// Set combustion and vapor pointers
template<class MixtureType>
void Foam::hsTwophaseMixtureThermo<MixtureType>::setPtrs
(
    combustionModels::rhoChemistryCombustionModel* combustion
)
{
    combustionPtr_ = combustion;
    
    //Set combustion pointer in all phases (to access reaction rates)
    liquid_.setCombustionPtr( combustion );
    vapor_.setCombustionPtr( combustion );
}


// Return field used for adaptive mesh refinement criteria
template<class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::hsTwophaseMixtureThermo<MixtureType>::getRefinementField() const
{
    // Refine ALL cells at the interface to the maximum level
    tmp<volScalarField> tRefinementField = interface_.getRefinementField();


    // Then enforce other, less-strict requirements in the rest of the domain
    /*forAll(this->Y(), i)
    {
        const volScalarField& Yi = this->Y()[i];
        tRefinementField().internalField() = max
        (
            tRefinementField().internalField(), 
            mag(fvc::grad(Yi)) * Foam::pow(mesh_.V(),1.0/3.0)
        );
    }*/

    //Include total volume source into refinement criteria
    /*tRefinementField().internalField() = max
    (
        tRefinementField().internalField(), 
        Sv_tot_.internalField() * 1e-6 //TODO Scale by cell volume mesh_.V()
    );*/

    //Include curl criteria from Popinet (Gerris), scaled by 0.5
    tRefinementField().internalField() = max
    (
        tRefinementField().internalField(), 
        Foam::mag(fvc::curl(U_)) * Foam::pow(mesh_.V(),1.0/3.0) * 0.5
    );

    return tRefinementField;
}



template<class MixtureType>
void Foam::hsTwophaseMixtureThermo<MixtureType>::calculateSurfaceTension()
{
    //Calculate the surface tension field
    dimensionedScalar one("one",dimLength,1.0);
    tmp<volScalarField> mask = pos(Foam::mag(interface_.kappa())*one - 1e-4)
                             + neg(vapor_ - 0.5);
    sigma_ = liquid_.sigma(T_, mask());
}

template<class MixtureType>
void Foam::hsTwophaseMixtureThermo<MixtureType>::solve
(
    volScalarField& rho,          // <- global density field
    const surfaceScalarField& phi // <- global volume flux field
)
{
    //Correct phases (correct liquid viscosity model)
    liquid_.correct();

    //Solve for reaction rates
    //Info<< "Solving combustion" << endl;
    //combustionPtr_->correct();

    //Solve for evaporation rates
    //Info<< "Solving evaporation" << endl;
    //calcEvaporation();

    //Do solving for phase volume fractions
    const Time& runTime = mesh_.time();
    const dictionary& pimpleDict = mesh_.solutionDict().subDict("PIMPLE");
    label nAlphaSubCycles(readLabel(pimpleDict.lookup("nAlphaSubCycles")));
    //scalar cAlpha(readScalar(pimpleDict.lookup("cAlpha")));
    //label nAlphaCorr(readLabel(pimpleDict.lookup("nAlphaCorr")));


    Info<< "Beginning alpha subcycle" << endl;
    if (nAlphaSubCycles > 1)
    {
        dimensionedScalar totalDeltaT = runTime.deltaT();
        
        // In subcycling mode mass fluxes are zeroed and incremented
        rhoPhi_ *= 0.0;
        liquid_.rhoPhi() *= 0.0;
        vapor_.rhoPhi() *= 0.0;
        
        for
        (
            subCycle<volScalarField> alphaSubCycle(liquid_, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            interface_.advect( liquid_.Sv_evap(), phi );
        
            //Update mass flux fields in each phase: phi = phiAlpha * rhof
            scalar fraction = (runTime.deltaT().value()/totalDeltaT.value());
            
            // Add the fractional mass fluxes to each phase
            liquid_.setRhoPhi( interface_.phiAlphaLiquid(), fraction );
            vapor_.setRhoPhi( interface_.phiAlphaVapor(), fraction );
            
            // total mass flux is the sum of the phase mass fluxes            
            rhoPhi_ += vapor_.rhoPhi() + liquid_.rhoPhi();
        }
    }
    else
    {
        interface_.advect( liquid_.Sv_evap(), phi );
            
        //Update mass flux fields in each phase: phi = phiAlpha * rhof
        liquid_.setRhoPhi( interface_.phiAlphaLiquid() );
        vapor_.setRhoPhi( interface_.phiAlphaVapor() );    
        
        // total mass flux is the sum of the phase mass fluxes
        rhoPhi_ = vapor_.rhoPhi() + liquid_.rhoPhi();  
    }

    // Update phase densities to satisfy continuity
    //   right now this just checks how far off they are, but does not change them
    liquid_.updateRho( interface_ );
    vapor_.updateRho( interface_ );
            
            
    
    
    //update global density field to satisfy continuity with new rhoPhi
    Foam::solve( fvm::ddt(rho) + fvc::div(rhoPhi_) );
    
    Info<< "Min, max rho = " << Foam::gMin(rho) << ", "
        << Foam::gMax(rho) << endl;
    
    
    
    
    //set interface and domain-related source terms
    //tmp<surfaceScalarField> tw = interface_.scTransferWeights();
    
    //no actual transfers to do yet
    
    
    
    

/*
    tmp<surfaceVectorField> ucL = liquid_.calculateDs
    (
        mesh_.time().time().value() > phaseRelaxTime_, 
        combustionPtr_->turbulence(),
        interface_.alphaLiquidf()
    );
    
    tmp<surfaceVectorField> ucV = vapor_.calculateDs
    (
        mesh_.time().time().value() > phaseRelaxTime_, 
        combustionPtr_->turbulence(),
        interface_.alphaVaporf()
    );
    
    Info<< "Max ucL = " << Foam::max(Foam::mag(ucL())).value() << endl;
    Info<< "Max ucV = " << Foam::max(Foam::mag(ucV())).value() << endl;

    // Create convection scheme for species
    tmp<fv::convectionScheme<scalar> > mvConvection
    (
        fv::convectionScheme<scalar>::New
        (
            mesh_,
            fields_,
            rhoPhi_,
            mesh_.divScheme("div(phi,Yi_h)")
        )
    );
    
    // Solve for subspecie transport within each phase
    scalar FoLiq = liquid_.solveSubSpecies(rho, rhoPhi_, p_, T_, liquid_, ucL(), mvConvection, interface_.alphaLiquidf());
    scalar FoVap = vapor_.solveSubSpecies(rho, rhoPhi_, p_, T_, liquid_, ucV(), mvConvection, interface_.alphaVaporf());

    scalar MaxFo = (FoVap > FoLiq) ? FoVap : FoLiq;
    
    // Calculate Ysum
    YSum_ = liquid_.Yp() + vapor_.Yp();
    
    Foam::Info << "Min,Max Ysum = " << Foam::min(YSum_).value() 
               << ", " << Foam::max(YSum_).value() << Foam::endl;
               
    return MaxFo;*/

}

template<class MixtureType>
Foam::tmp<Foam::volScalarField> 
Foam::hsTwophaseMixtureThermo<MixtureType>::S_evap() const
{
    return liquid_.S_evap(p_,T_);
}


template<class MixtureType>
void Foam::hsTwophaseMixtureThermo<MixtureType>::correct()
{
    if (debug)
    {
        Info<<"entering hsTwophaseMixtureThermo<MixtureType>::correct()"<<endl;
    }

    calculate();

    if (debug)
    {
        Info<<"exiting hsTwophaseMixtureThermo<MixtureType>::correct()"<<endl;
    }
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::hsTwophaseMixtureThermo<MixtureType>::hc() const
{
    tmp<volScalarField> thc
    (
        new volScalarField
        (
            IOobject
            (
                "hc",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            hs_.dimensions()
        )
    );

    volScalarField& hcf = thc();
    scalarField& hcCells = hcf.internalField();

    forAll(hcCells, celli)
    {
        hcCells[celli] = this->cellMixture(celli).Hc();
    }

    forAll(hcf.boundaryField(), patchi)
    {
        scalarField& hcp = hcf.boundaryField()[patchi];

        forAll(hcp, facei)
        {
            hcp[facei] = this->patchFaceMixture(patchi, facei).Hc();
        }
    }

    return thc;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::hsTwophaseMixtureThermo<MixtureType>::hs
(
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> ths(new scalarField(T.size()));
    scalarField& hs = ths();

    forAll(T, celli)
    {
        hs[celli] = this->cellMixture(cells[celli]).Hs(T[celli]);
    }

    return ths;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::hsTwophaseMixtureThermo<MixtureType>::hs
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> ths(new scalarField(T.size()));
    scalarField& hs = ths();

    forAll(T, facei)
    {
        hs[facei] = this->patchFaceMixture(patchi, facei).Hs(T[facei]);
    }

    return ths;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::hsTwophaseMixtureThermo<MixtureType>::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCp(new scalarField(T.size()));

    scalarField& cp = tCp();

    forAll(T, facei)
    {
        cp[facei] = this->patchFaceMixture(patchi, facei).Cp(T[facei]);
    }

    return tCp;
}



template<class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::hsTwophaseMixtureThermo<MixtureType>::Cp() const
{
    tmp<volScalarField> tCp
    (
        new volScalarField
        (
            IOobject
            (
                "Cp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimEnergy/dimMass/dimTemperature
        )
    );

    volScalarField& cp = tCp();

    scalarField& cpCells = cp.internalField();
    const scalarField& TCells = T_.internalField();

    forAll(TCells, celli)
    {
        cpCells[celli] = this->cellMixture(celli).Cp(TCells[celli]);
    }

    forAll(T_.boundaryField(), patchi)
    {
        cp.boundaryField()[patchi] = Cp(T_.boundaryField()[patchi], patchi);
    }

    return tCp;
}



template<class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::hsTwophaseMixtureThermo<MixtureType>::Cv
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCv(new scalarField(T.size()));

    scalarField& cv = tCv();

    forAll(T, facei)
    {
        cv[facei] = this->patchFaceMixture(patchi, facei).Cv(T[facei]);
    }

    return tCv;
}



template<class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::hsTwophaseMixtureThermo<MixtureType>::Cv() const
{
    tmp<volScalarField> tCv
    (
        new volScalarField
        (
            IOobject
            (
                "Cv",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimEnergy/dimMass/dimTemperature
        )
    );

    volScalarField& cv = tCv();

    scalarField& cvCells = cv.internalField();
    const scalarField& TCells = T_.internalField();

    forAll(TCells, celli)
    {
        cvCells[celli] = this->cellMixture(celli).Cv(TCells[celli]);
    }

    forAll(T_.boundaryField(), patchi)
    {
        cv.boundaryField()[patchi] = Cv(T_.boundaryField()[patchi], patchi);
    }

    return tCv;
}




template<class MixtureType>
bool Foam::hsTwophaseMixtureThermo<MixtureType>::read()
{
    if (hsReactionThermo::read())
    {
        MixtureType::read(*this);
        return true;
    }
    else
    {
        return false;
    }
}

template<class MixtureType> 
Foam::tmp<Foam::volScalarField> 
Foam::hsTwophaseMixtureThermo<MixtureType>::kappaV() const
{
    tmp<volScalarField> tkappa
    (
        new volScalarField
        (
            IOobject
            (
                "tkappa",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("tkappa", dimPower/dimLength/dimTemperature, 0.0)
        )
    );

    scalarField& kappaCells = tkappa().internalField();
    const scalarField& TCells = T_.internalField();

    forAll(TCells, celli)
    {
        kappaCells[celli] = this->cellMixture(celli).kappa( TCells[celli] );
    }

    forAll(T_.boundaryField(), patchi)
    {
        tkappa().boundaryField()[patchi] = kappa(T_.boundaryField()[patchi], patchi);
    }

    return tkappa;
}

template<class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::hsTwophaseMixtureThermo<MixtureType>::kappa
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tkappa(new scalarField(T.size()));

    scalarField& kappa = tkappa();

    forAll(T, facei)
    {
        kappa[facei] = this->patchFaceMixture(patchi, facei).kappa(T[facei]);
    }

    return tkappa;
}



template<class MixtureType>
tmp<volScalarField> 
Foam::hsTwophaseMixtureThermo<MixtureType>::kByCv
(
    const volScalarField& muEff
) const
{
   // kByCv = (alpha1*k1/Cv1 + alpha2*k2/Cv2) + muEff

    return liquid_*liquid_.k()/liquid_.Cv(T_)
         + vapor_*kappaV()/vapor_.Cv(T_)
         + muEff;
}

template<class MixtureType>
tmp<volScalarField> Foam::hsTwophaseMixtureThermo<MixtureType>::rCv() const
{
    // mixture.rCv() = (alphaL/CvL + alphaV/CvV)  
    return liquid_/liquid_.Cv(T_) + vapor_/vapor_.Cv(T_);
}



template<class MixtureType>
void Foam::hsTwophaseMixtureThermo<MixtureType>::setHs
(
    const volScalarField& T
)
{
    //Set hs based on an input T field
    // This is a backhanded way of setting T in an hsThermo object
    scalarField& hsCells = hs_.internalField();
    const scalarField& TCells = T.internalField();

    forAll(hsCells, celli)
    {
        hsCells[celli] = this->cellMixture(celli).Hs(TCells[celli]);
    }

    forAll(hs_.boundaryField(), patchi)
    {
        hs_.boundaryField()[patchi] == hs(T.boundaryField()[patchi], patchi);
    }

    hBoundaryCorrection(hs_);
}


// ************************************************************************* //
