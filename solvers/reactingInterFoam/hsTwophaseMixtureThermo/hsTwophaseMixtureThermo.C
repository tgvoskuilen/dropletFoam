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
    scalarField& alphaCells = alpha_.internalField();
    
    Info << "Calculating mixture temperature from hs" << endl;

    forAll(TCells, celli)
    {
        // references multiComponentMixture.C
        const typename MixtureType::thermoType& mixture =
            this->cellMixture(celli);

        TCells[celli] = mixture.THs(hsCells[celli], TCells[celli]);
        alphaCells[celli] = mixture.alpha(TCells[celli]);
    }

    forAll(T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pT = T_.boundaryField()[patchi];
        fvPatchScalarField& phs = hs_.boundaryField()[patchi];
        fvPatchScalarField& palpha = alpha_.boundaryField()[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture =
                    this->patchFaceMixture(patchi, facei);

                phs[facei] = mixture.Hs(pT[facei]);

                palpha[facei] = mixture.alpha(pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture =
                    this->patchFaceMixture(patchi, facei);

                pT[facei] = mixture.THs(phs[facei], pT[facei]);

                palpha[facei] = mixture.alpha(pT[facei]);
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
    alphaVapor_.correct(p_,T_);
    alphaLiquid_.correct(p_,T_);
    
    mu_ = alphaVapor_ * muV() + 
            alphaLiquid_ * alphaLiquid_.mu(p_, T_);
    
}



template<class MixtureType>
void Foam::hsTwophaseMixtureThermo<MixtureType>::calculatePsi()
{    
    psi_ = alphaVapor_*alphaVapor_.psi(T_);
    psi_.correctBoundaryConditions();
}


template<class MixtureType>
void Foam::hsTwophaseMixtureThermo<MixtureType>::calcEvaporation()
{
    Foam::Info << "Calculating evaporation zones" << Foam::endl;
    
    if( alphaLiquid_.subSpecies().size() > 1 )
    {
        tmp<volScalarField> tgradProd
        (
            new volScalarField
            (
                IOobject
                (
                    "tgradProd",
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar("tgradProd", dimless, 1.0)
            )
        );
        
        dimensionedScalar one("one",dimLength,1.0);
    
        forAllIter(PtrDictionary<subSpecie>, alphaLiquid_.subSpecies(), specieI)
        {
            tgradProd() *= Foam::mag( fvc::grad(specieI().Y())*one );
        }
        
        Info<< "Max grad prod = " << Foam::max(tgradProd()).value() << endl;
        
        evap_mask_ = pos(neg(tgradProd() - 100000.0)
                         + pos(alphaVaporSharp_ - 0.9) - SMALL);
    }
    else
    {
        evap_mask_ = pos(alphaVaporSharp_+1.0); //no multi-liquid mask
    }
    
    
    
    forAllIter(PtrDictionary<subSpecie>, alphaLiquid_.subSpecies(), specieI)
    {
        if (specieI().hasEvaporation())
        {
            specieI().evapModel().calculate
            (
                evap_mask_
            );
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
    alphaVapor_
    (
        "Vapor",
        subDict("Vapor"),
        mesh,
        this->Y(),
        this->speciesData()
    ),
    alphaLiquid_
    (
        "Liquid",
        subDict("Liquid"),
        mesh,
        this->Y(),
        this->speciesData()
    ),
    rhoPhi_
    (
        IOobject
        (
            "rho*phi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rho*phi", dimMass/dimTime, 0.0)
    ),
    phi_( rhoPhi_.db().lookupObject<surfaceScalarField>("phi") ),
    U_( rhoPhi_.db().lookupObject<volVectorField>("U") ),
    alphaVaporSharp_
    (
        IOobject
        (
            "alphaVaporSharp",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("alphaVaporSharp", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    kappaI_
    (
        IOobject
        (
            "kappaI",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("kappaI", dimless/dimLength, 0.0),
        zeroGradientFvPatchScalarField::typeName
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
    evap_mask_
    (
        IOobject
        (
            "evap_mask",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("evap_mask", dimless, 0.0),
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
    phaseRelaxTime_(readScalar(lookup("relaxTime")))
{

    // Create convection field table
    forAll(this->Y(), i)
    {
        fields_.add(this->Y()[i]);
    }

    //Set evaporation models
    Info<< "Setting evaporation models" << endl;

    forAllIter(PtrDictionary<subSpecie>, alphaLiquid_.subSpecies(), specieI)
    {
        if(specieI().hasEvaporation())
        {
            specieI().evapPtr() = evaporationModel::New
            (
                specieI().dict(),
                p_,
                T_,
                alphaLiquid_,
                alphaVapor_
            );
        }
    }

    calculateRho();
    rho_.oldTime();   
    
    rhoPhi_ = phi_ * fvc::interpolate(rho_);
    
    setHs(T_);

    calculate();
    
    alphaVapor_.setSpecies( rho_ );
    alphaLiquid_.setSpecies( rho_ );
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
    rho_ = alphaLiquid_*alphaLiquid_.rho(p_,T_) + 
            alphaVapor_*alphaVapor_.rho(p_,T_);
    rho_.correctBoundaryConditions();
    
    if( Foam::min(rho_).value() < SMALL )
    {
        Foam::Info<< "Min rho = " << Foam::min(rho_) << Foam::endl;
        Foam::Info<< "Min rho.if = " << Foam::min(rho_.internalField()) << Foam::endl;
        Foam::Info<< "Min p = " << Foam::min(p_) << Foam::endl;
        
        tmp<volScalarField> trhoV = alphaVapor_.rho(p_,T_);
        Foam::Info<< "Min rhoV = " << Foam::min(trhoV()) << Foam::endl;
        Foam::Info<< "Min rhoV.if = " << Foam::min(trhoV().internalField()) << Foam::endl;
        
        FatalErrorIn
        (
            "multiphaseReactingMixture::calculateRho()"
        )   << "Encountered a zero mixture density"
            << "\n    Check that phase boundary conditions have been applied"
            << "\n    and that a nonzero alpha is specified on each boundary."
            << exit(FatalError);
    }

    rho_.correctBoundaryConditions();
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

    forAllConstIter(PtrDictionary<subSpecie>, alphaVapor_.subSpecies(), specieI)
    {
        tden() += specieI().Y() / specieI().W();
    }
    
    forAllConstIter(PtrDictionary<subSpecie>, alphaLiquid_.subSpecies(), specieI)
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
        fvc::interpolate(alphaVapor_)*fvc::interpolate(muV())
      + fvc::interpolate(alphaLiquid_)*fvc::interpolate(alphaLiquid_.mu(p_,T_));

    return tmuf;
}


// Total enthalpy source/sink, including reactions and evaporation
template<class MixtureType> 
Foam::tmp<Foam::volScalarField> 
Foam::hsTwophaseMixtureThermo<MixtureType>::Sh() const
{
    return combustionPtr_->Sh() + alphaLiquid_.Sh_evap();
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


    tnearInt() = max(tnearInt(), pos(alphaLiquid_-lower)*pos(upper-alphaLiquid_));
    tnearInt() = max(tnearInt(), pos(alphaVapor_-lower)*pos(upper-alphaVapor_));

    return tnearInt;
}



// Face-centered surface tension force scalar
template<class MixtureType>
Foam::tmp<Foam::surfaceScalarField>
Foam::hsTwophaseMixtureThermo<MixtureType>::surfaceTensionForce()
{
    tmp<surfaceScalarField> tstf
    (
        new surfaceScalarField
        (
            IOobject
            (
                "surfaceTensionForce",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar
            (
                "surfaceTensionForce",
                dimensionSet(1, -2, -2, 0, 0),
                0.0
            )
        )
    );

    surfaceScalarField& stf = tstf();

    stf = fvc::interpolate(sigma_ * kappaI_) * fvc::snGrad(alphaVaporSharp_);

    return tstf;
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
    alphaLiquid_.setCombustionPtr( combustion );
    alphaVapor_.setCombustionPtr( combustion );
}


// Return field used for adaptive mesh refinement criteria
template<class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::hsTwophaseMixtureThermo<MixtureType>::getRefinementField() const
{
    tmp<volScalarField> tRefinementField
    (
        new volScalarField
        (
            IOobject
            (
                "tRefinementField",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("tRefinementField", dimless, 0.0)
        )
    );

    // Normalized Gradient Method
    //  RF = max( |grad alpha| * V^(1/3) )
    tRefinementField().internalField() = max
    (
        tRefinementField().internalField(), 
        20.0 * mag(fvc::grad(alphaLiquid_)) * Foam::pow(mesh_.V(),1.0/3.0)
    );


    forAll(this->Y(), i)
    {
        const volScalarField& Yi = this->Y()[i];
        tRefinementField().internalField() = max
        (
            tRefinementField().internalField(), 
            mag(fvc::grad(Yi)) * Foam::pow(mesh_.V(),1.0/3.0)
        );
    }

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
void Foam::hsTwophaseMixtureThermo<MixtureType>::calculateAlphaVapor()
{
    //Sharpen the remaining field
    scalar Cpc = 0.02;
    alphaVaporSharp_ = (Foam::min(Foam::max(alphaVapor_, 0.5*Cpc),1.0-0.5*Cpc)
                         - 0.5*Cpc)/(1.0-Cpc);
    alphaVaporSharp_.correctBoundaryConditions();
}




template<class MixtureType>
void Foam::hsTwophaseMixtureThermo<MixtureType>::calculateSurfaceTension()
{
    //Calculate interface curvature field
    // Cell gradient of alpha
    const volVectorField gradAlpha(fvc::grad(alphaVaporSharp_));

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

    // Face unit interface normal
    surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));

    // Simple expression for curvature
    kappaI_ = -fvc::div(nHatfv & mesh_.Sf());


    //Calculate the surface tension field
    dimensionedScalar one("one",dimLength,1.0);
    tmp<volScalarField> mask = pos(Foam::mag(kappaI_)*one - 1e-4)
                             + neg(alphaVaporSharp_ - 0.5);
    sigma_ = alphaLiquid_.sigma(T_, mask());
}

template<class MixtureType>
scalar Foam::hsTwophaseMixtureThermo<MixtureType>::solve
(
    volScalarField& rho
)
{
    calculateAlphaVapor();

    //Correct phases (correct liquid viscosity model)
    //alphaLiquid_.correct();

    //Solve for reaction rates
    Info<< "Solving combustion" << endl;
    combustionPtr_->correct();

    //Solve for evaporation rates
    Info<< "Solving evaporation" << endl;
    calcEvaporation();


    //Do solving for phase volume fractions
    const Time& runTime = mesh_.time();
    const dictionary& pimpleDict = mesh_.solutionDict().subDict("PIMPLE");
    label nAlphaSubCycles(readLabel(pimpleDict.lookup("nAlphaSubCycles")));
    scalar cAlpha(readScalar(pimpleDict.lookup("cAlpha")));
    label nAlphaCorr(readLabel(pimpleDict.lookup("nAlphaCorr")));


    Info<< "Beginning alpha subcycle" << endl;
    if (nAlphaSubCycles > 1)
    {
        surfaceScalarField rhoPhiSum(0.0*rhoPhi_);
        dimensionedScalar totalDeltaT = runTime.deltaT();

        for
        (
            subCycle<volScalarField> alphaSubCycle(alphaLiquid_, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            solveAlphas(cAlpha, nAlphaCorr);
            rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi_;
        }

        rhoPhi_ = rhoPhiSum;
    }
    else
    {
        solveAlphas(cAlpha, nAlphaCorr);
    }


    calculateAlphaVapor();
    calculateSurfaceTension();

    // Update density field to satisfy continuity with new mass flux field
    Foam::solve( fvm::ddt(rho) + fvc::div(rhoPhi_) );
    Info<< "Min,max rho = " << Foam::min(rho).value() << ", " 
        << Foam::max(rho).value() << endl;

    tmp<volVectorField> ucL = alphaLiquid_.calculateDs
    (
        0.95, 
        false, //mesh_.time().time().value() > phaseRelaxTime_, 
        combustionPtr_->turbulence()
    );
    
    tmp<volVectorField> ucV = alphaVapor_.calculateDs
    (
        0.95, 
        false, //mesh_.time().time().value() > phaseRelaxTime_, 
        combustionPtr_->turbulence()
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
    scalar FoLiq = alphaLiquid_.solveSubSpecies(rho, rhoPhi_, p_, T_, alphaLiquid_, ucL(), mvConvection);
    scalar FoVap = alphaVapor_.solveSubSpecies(rho, rhoPhi_, p_, T_, alphaLiquid_, ucV(), mvConvection);

    scalar MaxFo = (FoVap > FoLiq) ? FoVap : FoLiq;
    
    // Calculate Ysum
    YSum_ = alphaLiquid_.Yp() + alphaVapor_.Yp();
    
    Foam::Info << "Min,Max Ysum = " << Foam::min(YSum_).value() 
               << ", " << Foam::max(YSum_).value() << Foam::endl;
               
    return MaxFo;
}

template<class MixtureType>
Foam::tmp<Foam::volScalarField> 
Foam::hsTwophaseMixtureThermo<MixtureType>::S_evap() const
{
    return alphaLiquid_.S_evap(p_,T_);
}

template<class MixtureType>
void Foam::hsTwophaseMixtureThermo<MixtureType>::solveAlphas
(
    const scalar cAlpha,
    const label nAlphaCorr
)
{
    Info << "Starting solveAlphas" << endl;
    
    // For comparison with the code from compressibleInterFoam
    //  alpha1 = Liquid
    //  alpha2 = Vapor
    // however since we regard the liquid as incompressible considerable
    // simplifications can be made regarding the DpDt (dgdt) term
    
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");
    
    volScalarField divU(fvc::div(phi_));

    surfaceScalarField phic(mag(phi_/mesh_.magSf()));
    phic = min(cAlpha*phic, max(phic));
    
    //Calculate interface curvature field
    // Cell gradient of alpha
    const volVectorField gradAlpha(fvc::grad(alphaLiquid_)); //TODO: Use alphaVapor ?

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

    // Face unit interface normal
    surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));

    surfaceScalarField phir(phic*(nHatfv & mesh_.Sf()));
    
    for (int gCorr=0; gCorr<nAlphaCorr; gCorr++)
    {
        volScalarField::DimensionedInternalField Su
        (
            IOobject
            (
                "Su",
                mesh_.time().timeName(),
                mesh_
            ),
            // Divergence term is handled explicitly to be
            // consistent with the explicit transport solution
            divU*min(alphaLiquid_, scalar(1)) - alphaLiquid_.Sv_evap()
        );        

        surfaceScalarField phiAlphaL
        (
            fvc::flux
            (
                phi_,
                alphaLiquid_,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, alphaVapor_, alpharScheme),
                alphaLiquid_,
                alpharScheme
            )
        );

        MULES::explicitSolve
        (
            geometricOneField(),
            alphaLiquid_,
            phi_,
            phiAlphaL,
            zeroField(),
            Su,
            1,
            0
        );
        
        alphaLiquid_.max(0.0);
        
        // WARNING:
        //   This solution provides a rho and rhoPhi pair that do not satisfy
        //   continuity (DaDt is solved). The rhoEqn must be solved to get
        //   an appropriate density field before other equations in
        //   conservative form can be solved. Updating rho using
        //   rho = sum(alphaI * rhoI) will NOT work properly.

        surfaceScalarField rhoVf(fvc::interpolate(alphaVapor_.rho(p_,T_)));
        surfaceScalarField rhoLf(fvc::interpolate(alphaLiquid_.rho(p_,T_)));
        
        alphaLiquid_.rhoPhiAlpha() = phiAlphaL * rhoLf;
        alphaVapor_.rhoPhiAlpha() = phi_*rhoVf - phiAlphaL * rhoVf;

        rhoPhi_ = phiAlphaL*(rhoLf - rhoVf) + phi_*rhoVf;

        alphaVapor_ == scalar(1) - alphaLiquid_;
        alphaVapor_.max(0.0);
    }
    
    
    
    Info<< "Liquid phase volume fraction = "
        << alphaLiquid_.weightedAverage(mesh_.V()).value()
        << "  Min,max alphaV = " << min(alphaVapor_).value() 
        << ", " << max(alphaVapor_).value()
        << "  Min,max alphaL = " << min(alphaLiquid_).value()
        << ", " << max(alphaLiquid_).value()
        << endl;
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

    return alphaLiquid_*alphaLiquid_.k()/alphaLiquid_.Cv(T_)
         + alphaVapor_*kappaV()/alphaVapor_.Cv(T_)
         + muEff;
}

template<class MixtureType>
tmp<volScalarField> Foam::hsTwophaseMixtureThermo<MixtureType>::rCv() const
{
    // mixture.rCv() = (alphaL/CvL + alphaV/CvV)  
    return alphaLiquid_/alphaLiquid_.Cv(T_) + alphaVapor_/alphaVapor_.Cv(T_);
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
