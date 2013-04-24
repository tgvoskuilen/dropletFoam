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
    
    // calculates area_ from alphaL
    correctInterface();
    
    // uses area_ and alpha to set mask
    alphaLiquid_.setPhaseMasks(phaseMaskTol_, p_, T_, area_ );
    alphaVapor_.setPhaseMasks(phaseMaskTol_, p_, T_, area_ );
    
    // uses mask to set rhoAlpha
    alphaVapor_.correct(p_,T_);
    alphaLiquid_.correct(p_,T_);
    
    mu_ = alphaVapor_.sharp(0.0)*muV() + alphaLiquid_.sharp(0.0)*alphaLiquid_.mu(p_, T_);
    mu_.correctBoundaryConditions();
    
    psi_ = alphaVapor_.sharp(0.0)*alphaVapor_.psi(T_);
    psi_.correctBoundaryConditions();
    
    rho_ = alphaLiquid_.sharp(0.0)*alphaLiquid_.rho(p_,T_) + 
            alphaVapor_.sharp(0.0)*alphaVapor_.rho(p_,T_);
    rho_.correctBoundaryConditions();
}

template<class MixtureType>
void Foam::hsTwophaseMixtureThermo<MixtureType>::correctInterface()
{
    //Update kappa, sigma_, alphaSmooth, area_
    
    
    //Smooth alphavapor
    alphaVaporSmooth_ = alphaVapor_.sharp(0.1);
    dimensionedScalar dA = pow(min(mesh_.V()),2.0/3.0)/30.0;
    
    for( label i = 0; i < 15; ++i )
    {
        alphaVaporSmooth_ += dA * fvc::laplacian(alphaVaporSmooth_);
    }
    alphaVaporSmooth_.min(1.0);
    alphaVaporSmooth_.max(0.0);
    alphaVaporSmooth_.correctBoundaryConditions();
    
    
    /*
    //Sharpen the remaining field
    scalar Cpc = 0.02;
    alphaVaporSharp_ = (Foam::min(Foam::max(alphaVapor_, 0.5*Cpc),1.0-0.5*Cpc)
                         - 0.5*Cpc)/(1.0-Cpc);
    alphaVaporSharp_.correctBoundaryConditions();
    */
    
    
    
    //Calculate interface curvature field
    // Cell gradient of alphaVaporSmooth
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
    
    

    
    
    //Calculate area
    //hardt method
    tmp<volScalarField> C = alphaLiquid_.sharp(0.01);
    
    dimensionedScalar eps("eps",dimArea,SMALL);
    tmp<volScalarField> Cp = Foam::mag(fvc::grad(C()));
    dimensionedScalar N = fvc::domainIntegrate(Cp()) 
        / (fvc::domainIntegrate((1-C())*(1-C())*Cp()) + eps);

    area_ = N*(1-C())*(1-C())*Cp;
    //area_ = Cp;
    
    const volScalarField::DimensionedInternalField& V = mesh_.V();
    area_.dimensionedInternalField() *= pos
    (
        area_.dimensionedInternalField()
      - Foam::pow(V,-1.0/3.0)/150.0
    );
    area_.correctBoundaryConditions();
    
    
}

template<class MixtureType>
void Foam::hsTwophaseMixtureThermo<MixtureType>::calcPhaseChange()
{
    Foam::Info << "Calculating phase change zones" << Foam::endl;
    
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
        
        // only allow evaporation if gradprod < 1e5 OR alphav > 0.5
        evap_mask_ = pos(neg(tgradProd() - 100000.0)
                         + pos(alphaVapor_ - 0.5) - SMALL);
    }
    else
    {
        evap_mask_ = pos(alphaVapor_+1.0); //no multi-liquid mask
    }
    
    
    forAllIter(PtrDictionary<mixturePhaseChangeModel>, phaseChangeModels_, pcmI)
    {
        pcmI().calculate
        (
            evap_mask_, //legacy, maybe remove
            area_
        );
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
    phaseChangeModels_
    (
        lookup("phaseChangeReactions"),
        mixturePhaseChangeModel::iNew
        (
            mesh,
            alphaLiquid_,
            alphaVapor_,
            this->speciesData()
        )
    ),
    area_
    (
        IOobject
        (
            "phaseArea",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("area", dimless/dimLength, 0.0)
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
    divPhaseChange_
    (
        IOobject
        (
            "divPhaseChange",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("divPhaseChange", dimless/dimTime, 0.0)
    ),
    divComp_
    (
        IOobject
        (
            "divComp",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::DDt(phi_,p_)/(p_+dimensionedScalar("ps",dimPressure,1))
    ),
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
    divPhaseChange_.oldTime();
    divComp_.oldTime();
    alphaLiquid_.setOtherPhase( &alphaVapor_ );
    alphaVapor_.setOtherPhase( &alphaLiquid_ );

    setHs();
    calculate();
    
    rhoPhi_ = phi_ * fvc::interpolate(rho_);
    rho_.oldTime();
    
    //Set Yp values and allow slight diffusion to expand to phase mask boundary
    alphaVapor_.setSpecies( alphaLiquid_.rhoAlpha() );
    alphaLiquid_.setSpecies( alphaVapor_.rhoAlpha() );
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


template<class MixtureType> 
Foam::tmp<Foam::volScalarField> 
Foam::hsTwophaseMixtureThermo<MixtureType>::dQ_phaseChange() const
{
    tmp<volScalarField> dQ
    (
        new volScalarField
        (
            IOobject
            (
                "dQ",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("dQ", dimPower/dimVolume, 0.0)
        )
    );
    
    forAllConstIter
    (
        PtrDictionary<mixturePhaseChangeModel>, 
        phaseChangeModels_, 
        pcmI
    )
    {
        dQ() += pcmI().Sh();
    }
            
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

    //stf = fvc::interpolate(sigma_*kappaI_ + alphaVapor_.pRecoil()) * fvc::snGrad(alphaVapor_.H());
    
    //dimensionedScalar meanRho("meanRho",dimDensity,500); //TODO: calc mean rho
    
    //surfaceScalarField rhof = fvc::interpolate(rho_)/meanRho;
    
    stf = fvc::interpolate(sigma_*kappaI_) * fvc::snGrad(alphaVapor_.sharp(0.01));// * rhof;
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
        
    forAllIter
    (
        PtrDictionary<mixturePhaseChangeModel>, 
        phaseChangeModels_, 
        pcmI
    )
    {
        pcmI().setPtr( combustion );
    }
    
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
        1000.0 * mag(fvc::grad(alphaLiquid_.sharp(0.01))) * Foam::pow(mesh_.V(),1.0/3.0)
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
scalar Foam::hsTwophaseMixtureThermo<MixtureType>::solve
(
    volScalarField& rho
)
{
    //Solve for reaction rates
    Info<< "Solving combustion" << endl;
    rho_ = alphaLiquid_.rhoAlpha() + alphaVapor_.rhoAlpha();
    alphaLiquid_.updateGlobalYs( alphaVapor_.rhoAlpha() );
    alphaVapor_.updateGlobalYs( alphaLiquid_.rhoAlpha() );
    combustionPtr_->correct();

    //Solve for evaporation rates
    Info<< "Solving phase change" << endl;
    calcPhaseChange();

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
    
    // correct() updates rho_, psi_, mu_ and calls correct() on phases
    // Update properties and set rhoAlpha to satisfy continuity with the
    // new rhoPhiAlpha mass flux field
    correct();
    
    // Because continuity is satisfied on a per-phase basis, it should also
    // be satisfied on a global basis with rho and rhoPhi_
    rho = alphaLiquid_.rhoAlpha() + alphaVapor_.rhoAlpha();
    rhoPhi_ = alphaLiquid_.rhoPhiAlpha() + alphaVapor_.rhoPhiAlpha();
    
    Info<< "Min,max rho = " << Foam::min(rho_).value() << ", " 
        << Foam::max(rho_).value() << endl;
        
        
    YSum_ = alphaLiquid_.Yp() + alphaVapor_.Yp();
    
    Info<< "Min,Max Ysum = " << Foam::min(YSum_).value() 
        << ", " << Foam::max(YSum_).value() << endl;
    
    // Calculate diffusion coefficient for each phase
    alphaLiquid_.calculateDs( combustionPtr_->turbulence().muEff() );
    alphaVapor_.calculateDs( combustionPtr_->turbulence().muEff() );
       
    // Solve for subspecie transport within each phase
    scalar FoLiq = alphaLiquid_.solveSubSpecies( p_, T_, phaseChangeModels_ );
    scalar FoVap = alphaVapor_.solveSubSpecies( p_, T_, phaseChangeModels_ );

    // Update global mass fractions based on phase-based mass fractions
    alphaLiquid_.updateGlobalYs( alphaVapor_.rhoAlpha() );
    alphaVapor_.updateGlobalYs( alphaLiquid_.rhoAlpha() );
    
    //scalar MaxFo = (FoVap > FoLiq) ? FoVap : FoLiq;
    
    // Calculate Ysum
    /*YSum_ = alphaLiquid_.Yp() + alphaVapor_.Yp();
    PtrList<volScalarField>& Ys = this->Y();
    Info<< "Min,Max Ysum = " << Foam::min(YSum_).value() 
        << ", " << Foam::max(YSum_).value() << endl;
        
    forAll(Ys, i)
    {   
        Ys[i] /= YSum_;
    }*/
    
    YSum_ = alphaLiquid_.Yp() + alphaVapor_.Yp();
    
    Info<< "Min,Max Ysum = " << Foam::min(YSum_).value() 
        << ", " << Foam::max(YSum_).value() << endl;
    
    return 0.0; //MaxFo; //TODO: Use DiNum from cht
}



template<class MixtureType>
Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::hsTwophaseMixtureThermo<MixtureType>::TSuSp() const
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
                dimensionedScalar("TSp",dimPower/dimVolume/dimTemperature,0.0)
            )
        )
    );
    
    forAllConstIter
    (
        PtrDictionary<mixturePhaseChangeModel>, 
        phaseChangeModels_, 
        pcmI
    )
    {
        Pair<tmp<volScalarField> > pcmTSuSp = pcmI().TSuSp();
        tTSuSp.first()() += pcmTSuSp.first();
        tTSuSp.second()() += pcmTSuSp.second();
    }
    
    tTSuSp.first()() *= rCv();
    tTSuSp.second()() *= rCv();
    
    return tTSuSp;
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

    surfaceScalarField phir(phic*(nHatfv & mesh_.Sf()));
        
    surfaceScalarField phiAlphaL("phiAlphaL",phi_);
    
    divPhaseChange_ *= 0.0;
    volScalarField Sv
    (
        IOobject
        (
            "Sv",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("Sv",dimless/dimTime,0.0)
    );
    
    forAllConstIter
    (
        PtrDictionary<mixturePhaseChangeModel>, 
        phaseChangeModels_, 
        pcmI
    )
    {
        divPhaseChange_ += pcmI().Vdot("Vapor") + pcmI().Vdot("Liquid");
        Sv += pcmI().Vdot("Liquid");
    }
    
    
    for (int gCorr=0; gCorr<nAlphaCorr; gCorr++)
    {
        volScalarField Sp
        (
            IOobject
            (
                "Sp",
                mesh_.time().timeName(),
                mesh_
            ),
            -divPhaseChange_
        );
        
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
            divU*min(alphaLiquid_, scalar(1)) + Sv
        );
        
        forAll(divComp_, celli)
        {
            if (divComp_[celli] > 0.0 && alphaLiquid_[celli] > 0.0)
            {
                Sp[celli] -= divComp_[celli]*alphaLiquid_[celli];
                Su[celli] += divComp_[celli]*alphaLiquid_[celli];
            }
            else if (divComp_[celli] < 0.0 && alphaLiquid_[celli] < 1.0)
            {
                Sp[celli] += divComp_[celli]*(1.0 - alphaLiquid_[celli]);
            }
        }
                
        phiAlphaL = 
        (
            fvc::flux
            (
                phi_,
                alphaLiquid_,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, 1.0 - alphaLiquid_, alpharScheme),
                alphaLiquid_,
                alpharScheme
            )
        );


        //MULES::explicitSolve
        MULES::implicitSolve
        (
            geometricOneField(),
            alphaLiquid_,
            phi_,
            phiAlphaL,
            Sp,
            Su,
            1,
            0
        );
        
    }
    
    Info<< "Liquid phase volume fraction = "
        << alphaLiquid_.weightedAverage(mesh_.V()).value()
        << "  Min,max alphaV = " << min(alphaVapor_).value() 
        << ", " << max(alphaVapor_).value()
        << "  Min,max alphaL = " << min(alphaLiquid_).value()
        << ", " << max(alphaLiquid_).value()
        << endl;
        
    alphaLiquid_.max(0.0);
    //alphaLiquid_.min(1.0);
    
    surfaceScalarField rhoVf(fvc::interpolate(alphaVapor_.rho(p_,T_)));
    surfaceScalarField rhoLf(fvc::interpolate(alphaLiquid_.rho(p_,T_)));
    
    alphaLiquid_.rhoPhiAlpha() += f * phiAlphaL * rhoLf;
    alphaVapor_.rhoPhiAlpha() += f * (phi_ - phiAlphaL) * rhoVf;

    rhoPhi_ += f * (phiAlphaL*(rhoLf - rhoVf) + phi_*rhoVf);

    alphaVapor_ == scalar(1) - alphaLiquid_;
    alphaVapor_.max(0.0);
        
    // Re-sharpen alpha field
    //  This is not mass conserving, but prevents excessive floatsom
    //volScalarField& alphaL = alphaLiquid_;
    //alphaL = alphaLiquid_.sharp(1e-3);
    //alphaVapor_ == scalar(1) - alphaLiquid_;
        
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
            dimensionedScalar("tkappa", dimPower/dimLength/dimTemperature, 0.2)
        )
    );
/*
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
*/
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
