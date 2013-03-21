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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * *  //

template<class MixtureType>
void Foam::hsTwophaseMixtureThermo<MixtureType>::calculate()
{
    // Update the following stored and inherited fields:
    //   hs_, rho_, psi_, mu_, alpha_
    // These are functions of the stored p_ and T_ fields
    //
    const scalarField& hsCells = hs_.internalField();

    scalarField& TCells = T_.internalField();
    scalarField& alphaCells = alpha_.internalField();
    
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
    
    alphaVapor_.correct(p_,T_);
    alphaLiquid_.correct(p_,T_);
    
    mu_ = alphaVapor_*muV() + alphaLiquid_*alphaLiquid_.mu(p_, T_);
    mu_.correctBoundaryConditions();
    
    psi_ = alphaVapor_*alphaVapor_.psi(T_);
    psi_.correctBoundaryConditions();
    
    rho_ = alphaLiquid_*alphaLiquid_.rho(p_,T_) + 
            alphaVapor_*alphaVapor_.rho(p_,T_);
    rho_.correctBoundaryConditions();
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
    
        forAllIter
        (
            PtrDictionary<subSpecie>, alphaLiquid_.subSpecies(), specieI
        )
        {
            tgradProd() *= Foam::mag( fvc::grad(specieI().Y())*one );
        }
        
        Info<< "Max grad prod = " << Foam::max(tgradProd()).value() << endl;
        
        evap_mask_ = pos(neg(tgradProd() - 100000.0)
                         + pos(alphaVapor_ - 0.9) - SMALL);
    }
    else
    {
        evap_mask_ = pos(alphaVapor_+1.0); //no multi-liquid mask
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * *  //

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
    alphaVaporSmooth_
    (
        IOobject
        (
            "alphaVaporSmooth",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("alphaVaporSmooth", dimless, 0.0),
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
    phaseMaskTol_(readScalar(lookup("phaseMaskTol")))
{
    
    alphaLiquid_.setOtherPhase( &alphaVapor_ );
    alphaVapor_.setOtherPhase( &alphaLiquid_ );

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

    setHs();
    calculate();
    
    rhoPhi_ = phi_ * fvc::interpolate(rho_);
    rho_.oldTime();
    
    // Define phase boundary masks to prevent artificial cross-phase mixing
    //  Masks are based on current values of alphas
    alphaLiquid_.setPhaseMasks(phaseMaskTol_);
    alphaVapor_.setPhaseMasks(phaseMaskTol_);
    
    //Set Yp values and allow slight diffusion to expand to phase mask boundary
    alphaVapor_.setSpecies( rho_ );
    alphaLiquid_.setSpecies( rho_ );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class MixtureType>
Foam::hsTwophaseMixtureThermo<MixtureType>::~hsTwophaseMixtureThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //


// Cell centered mixture viscosity of the gas, calculated by Sutherland
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
        //TODO: Cell mixture of only gas species here, or apply mu mixture rule
        muCells[celli] = this->cellMixture(celli).mu( TCells[celli] );
    }

    forAll(T_.boundaryField(), patchi)
    {
        tmu().boundaryField()[patchi] = mu
        (
            T_.boundaryField()[patchi], 
            patchi
        );
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


// Total enthalpy source/sink, including reactions and evaporation
template<class MixtureType> 
Foam::tmp<Foam::volScalarField> 
Foam::hsTwophaseMixtureThermo<MixtureType>::Sh() const
{
    return combustionPtr_->Sh() + alphaLiquid_.Sh_evap();
}

template<class MixtureType> 
Foam::tmp<Foam::volScalarField> 
Foam::hsTwophaseMixtureThermo<MixtureType>::dQ_evap() const
{    
    tmp<volScalarField> dQ = alphaLiquid_.Sh_evap();

    dQ().dimensionedInternalField() *= mesh_.V();
    dQ().correctBoundaryConditions();

    return dQ;
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


    tnearInt() = max
    (
        tnearInt(), 
        pos(alphaLiquid_-lower)*pos(upper-alphaLiquid_)
    );

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

    stf = fvc::interpolate(sigma_ * kappaI_) * fvc::snGrad(alphaVapor_);

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
        1000.0 * mag(fvc::grad(alphaLiquid_)) * Foam::pow(mesh_.V(),1.0/3.0)
    );

/*
    forAll(this->Y(), i)
    {
        const volScalarField& Yi = this->Y()[i];
        tRefinementField().internalField() = max
        (
            tRefinementField().internalField(), 
            mag(fvc::grad(Yi)) * Foam::pow(mesh_.V(),1.0/3.0)
        );
    }
*/
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
    //Smooth alphavapor
    alphaVaporSmooth_ = alphaVapor_;
    dimensionedScalar dA = pow(min(mesh_.V()),2.0/3.0)/30.0;
    
    for( label i = 0; i < 15; ++i )
    {
        alphaVaporSmooth_ += dA * fvc::laplacian(alphaVaporSmooth_);
    }
    
    alphaVaporSmooth_.correctBoundaryConditions();
    
    
    /*
    //Sharpen the remaining field
    scalar Cpc = 0.02;
    alphaVaporSharp_ = (Foam::min(Foam::max(alphaVapor_, 0.5*Cpc),1.0-0.5*Cpc)
                         - 0.5*Cpc)/(1.0-Cpc);
    alphaVaporSharp_.correctBoundaryConditions();
    */
}




template<class MixtureType>
void Foam::hsTwophaseMixtureThermo<MixtureType>::calculateSurfaceTension()
{
    //Calculate interface curvature field
    // Cell gradient of alpha
    const volVectorField gradAlpha(fvc::grad(alphaVaporSmooth_));

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

    // Face unit interface normal
    surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));

    // Simple expression for curvature
    kappaI_ = -fvc::div(nHatfv & mesh_.Sf());


    //Calculate the surface tension field
    dimensionedScalar one("one",dimLength,1.0);
    tmp<volScalarField> mask = pos(Foam::mag(kappaI_)*one - 1e-4)
                             + neg(alphaVapor_ - 0.1);
    sigma_ = alphaLiquid_.sigma(T_, mask());
}



template<class MixtureType>
scalar Foam::hsTwophaseMixtureThermo<MixtureType>::solve
(
    volScalarField& rho
)
{
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

    // Zero out mass flux fields and increment them at each iteration
    rhoPhi_ *= 0.0;
    alphaLiquid_.rhoPhiAlpha() *= 0.0;
    alphaVapor_.rhoPhiAlpha() *= 0.0;
        
    Info<< "Beginning alpha subcycle" << endl;
    if (nAlphaSubCycles > 1)
    {
        scalar totalDeltaT = runTime.deltaTValue();
    
        for
        (
            subCycle<volScalarField> 
                alphaSubCycle(alphaLiquid_, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            scalar f = runTime.deltaTValue() / totalDeltaT;
            solveAlphas(cAlpha, nAlphaCorr, f);
        }
    }
    else
    {
        solveAlphas(cAlpha, nAlphaCorr);
    }

    // Update interface
    calculateAlphaVapor();
    calculateSurfaceTension();

    // Update density field to satisfy continuity with new mass flux field
    Foam::solve( fvm::ddt(rho) + fvc::div(rhoPhi_) );
    
    Info<< "Min,max rho = " << Foam::min(rho).value() << ", " 
        << Foam::max(rho).value() << endl;

    // Define phase boundary masks to prevent artificial cross-phase mixing
    alphaLiquid_.setPhaseMasks(phaseMaskTol_);
    alphaVapor_.setPhaseMasks(phaseMaskTol_);

    // Calculate diffusion coefficient for each phase
    alphaLiquid_.calculateDs( combustionPtr_->turbulence().muEff() );
    alphaVapor_.calculateDs( combustionPtr_->turbulence().muEff() );
       
    // Solve for subspecie transport within each phase
    scalar FoLiq = alphaLiquid_.solveSubSpecies( p_, T_ );
    scalar FoVap = alphaVapor_.solveSubSpecies( p_, T_ );

    // Update global mass fractions based on phase-based mass fractions
    alphaLiquid_.updateGlobalYs( alphaVapor_.rhoAlpha() );
    alphaVapor_.updateGlobalYs( alphaLiquid_.rhoAlpha() );
    
    scalar MaxFo = (FoVap > FoLiq) ? FoVap : FoLiq;
    
    // Calculate Ysum
    YSum_ = alphaLiquid_.Yp() + alphaVapor_.Yp();
    
    Info<< "Min,Max Ysum = " << Foam::min(YSum_).value() 
        << ", " << Foam::max(YSum_).value() << endl;
    
    return 0.0; //MaxFo; //TODO: Use DiNum from cht
}

template<class MixtureType>
Foam::tmp<Foam::volScalarField> 
Foam::hsTwophaseMixtureThermo<MixtureType>::S_evap() const
{
    return alphaLiquid_.S_evap(p_,T_);
}


template<class MixtureType>
Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::hsTwophaseMixtureThermo<MixtureType>::SuSp_evap() const
{
    return alphaLiquid_.pSuSp(p_,T_);
}



template<class MixtureType>
void Foam::hsTwophaseMixtureThermo<MixtureType>::solveAlphas
(
    const scalar cAlpha,
    const label nAlphaCorr,
    scalar f
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
    const volVectorField gradAlpha(fvc::grad(alphaLiquid_));

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

    // Face unit interface normal
    surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));

    //surfaceScalarField phiRecoil(fvc::interpolate(alphaVapor_.URecoil()) & mesh_.Sf());

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
        
        alphaLiquid_.rhoPhiAlpha() += f * phiAlphaL * rhoLf;
        alphaVapor_.rhoPhiAlpha() += f * (phi_ - phiAlphaL) * rhoVf;

        rhoPhi_ += f * (phiAlphaL*(rhoLf - rhoVf) + phi_*rhoVf);

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
        tkappa().boundaryField()[patchi] = 
            kappa(T_.boundaryField()[patchi], patchi);
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
void Foam::hsTwophaseMixtureThermo<MixtureType>::setHs()
{
    //Set hs based on an input T field
    // This is a backhanded way of setting T in an hsThermo object
    scalarField& hsCells = hs_.internalField();
    const scalarField& TCells = T_.internalField();

    forAll(hsCells, celli)
    {
        hsCells[celli] = this->cellMixture(celli).Hs(TCells[celli]);
    }

    forAll(hs_.boundaryField(), patchi)
    {
        hs_.boundaryField()[patchi] == hs(T_.boundaryField()[patchi], patchi);
    }

    hBoundaryCorrection(hs_);
}


// ************************************************************************* //
