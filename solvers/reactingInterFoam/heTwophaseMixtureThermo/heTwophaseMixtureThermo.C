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

#include "heTwophaseMixtureThermo.H"
#include "fvMesh.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * *  //

template<class ThermoType, class MixtureType>
void Foam::heTwophaseMixtureThermo<ThermoType,MixtureType>::calculate()
{
    // Update the following stored and inherited fields:
    //   he_, rho_, psi_, mu_, alpha_
    // These are functions of the stored p_ and T_ fields
    //
    const scalarField& hCells = this->he_.internalField();
    const scalarField& pCells = this->p_.internalField();
    
    scalarField& TCells = this->T_.internalField();
    //psi
    //rho
    //mu
    scalarField& alphaCells = this->alpha_.internalField();
    
    forAll(TCells, celli)
    {
        // references multiComponentMixture.C
        const typename MixtureType::thermoType& mixture =
            this->cellMixture(celli);
            
        TCells[celli] = mixture.THE
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );
        
        //psi
        //rho
        
        //mu
        alphaCells[celli] = mixture.alphah(pCells[celli], TCells[celli]);
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        //psi
        //rho
        
        fvPatchScalarField& ph = this->he_.boundaryField()[patchi];
        
        //mu
        fvPatchScalarField& palpha = this->alpha_.boundaryField()[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture =
                    this->patchFaceMixture(patchi, facei);

                ph[facei] = mixture.HE(pp[facei], pT[facei]);
                
                //psi
                //rho
                //mu
                palpha[facei] = mixture.alphah(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture =
                    this->patchFaceMixture(patchi, facei);

                pT[facei] = mixture.THE(ph[facei], pp[facei], pT[facei]);
                
                //psi
                //rho
                //mu
                palpha[facei] = mixture.alphah(pp[facei], pT[facei]);
            }
        }
    }
    
    // calculates area_ from alphaL
    correctInterface();
    
    // uses area_ and alpha to set mask
    alphaLiquid_.setPhaseMasks(phaseMaskTol_, this->p_, this->T_, phaseChangeModels_, area_ );
    alphaVapor_.setPhaseMasks(phaseMaskTol_, this->p_, this->T_, phaseChangeModels_, area_ );
    
    // uses mask to set rhoAlpha
    alphaVapor_.correct(this->p_,this->T_);
    alphaLiquid_.correct(this->p_,this->T_);
    
    
    //Now do psi, rho, and mu skipped above
    
    this->mu_ = alphaVapor_.sharp(0.0)*muV() + alphaLiquid_.sharp(0.0)*alphaLiquid_.mu(this->p_, this->T_);
    this->mu_.correctBoundaryConditions();
    
    this->psi_ = alphaVapor_.sharp(0.0)*alphaVapor_.psi(this->T_);
    this->psi_.correctBoundaryConditions();
    
    this->rho_ = alphaLiquid_.sharp(0.0)*alphaLiquid_.rho(this->p_,this->T_) + 
            alphaVapor_.sharp(0.0)*alphaVapor_.rho(this->p_,this->T_);
    this->rho_.correctBoundaryConditions();
}

template<class ThermoType, class MixtureType>
void Foam::heTwophaseMixtureThermo<ThermoType,MixtureType>::correctInterface()
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
    sigma_ = alphaLiquid_.sigma(this->T_, mask());
    
    

    
    
    //Calculate area
    //hardt method
    tmp<volScalarField> C = alphaLiquid_.sharp(0.0,0.01);
    
    word gradScheme("grad(alphaVaporSmooth)");
    
    dimensionedScalar eps("eps",dimArea,SMALL);
    tmp<volScalarField> Cp = Foam::mag(fvc::grad(C(),gradScheme));
    dimensionedScalar N = fvc::domainIntegrate(Cp()) 
        / (fvc::domainIntegrate((1-C())*(1-C())*Cp()) + eps);

    area_ = N*(1-C())*(1-C())*Cp;
    //area_ = Cp;
    
    /*const volScalarField::DimensionedInternalField& V = mesh_.V();
    area_.dimensionedInternalField() *= pos
    (
        area_.dimensionedInternalField()
      - Foam::pow(V,-1.0/3.0)/150.0
    );
    area_.correctBoundaryConditions();*/
    
    
}

template<class ThermoType, class MixtureType>
void Foam::heTwophaseMixtureThermo<ThermoType,MixtureType>::calcPhaseChange()
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

template<class ThermoType, class MixtureType>
Foam::heTwophaseMixtureThermo<ThermoType,MixtureType>::heTwophaseMixtureThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    heThermo<ThermoType,MixtureType>(mesh, phaseName),
    mesh_(mesh),
    combustionPtr_(NULL),
    alphaVapor_
    (
        "Vapor",
        ThermoType::subDict("Vapor"),
        mesh,
        this->Y(),
        this->speciesData()
    ),
    alphaLiquid_
    (
        "Liquid",
        ThermoType::subDict("Liquid"),
        mesh,
        this->Y(),
        this->speciesData()
    ),
    phaseChangeModels_
    (
        ThermoType::lookup("phaseChangeReactions"),
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
        fvc::DDt(phi_,this->p_)/(this->p_+dimensionedScalar("ps",dimPressure,1))
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
    phaseMaskTol_(readScalar(ThermoType::lookup("phaseMaskTol")))
{
    this->T_.correctBoundaryConditions();
    this->p_.correctBoundaryConditions();
    
    divPhaseChange_.oldTime();
    divComp_.oldTime();
    alphaLiquid_.setOtherPhase( &alphaVapor_ );
    alphaVapor_.setOtherPhase( &alphaLiquid_ );

    setHE(); //non-private copy of 'init' from heThermo
    calculate();
    
    rhoPhi_ = phi_ * fvc::interpolate(this->rho_);
    this->rho_.oldTime();
    
    //Set Yp values and allow slight diffusion to expand to phase mask boundary
    alphaVapor_.setSpecies( alphaLiquid_.rhoAlpha() );
    alphaLiquid_.setSpecies( alphaVapor_.rhoAlpha() );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class ThermoType, class MixtureType>
Foam::heTwophaseMixtureThermo<ThermoType,MixtureType>::~heTwophaseMixtureThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //


// Cell centered mixture viscosity of the gas, calculated by Sutherland
template<class ThermoType, class MixtureType> 
Foam::tmp<Foam::volScalarField> 
Foam::heTwophaseMixtureThermo<ThermoType,MixtureType>::muV() const
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
    const scalarField& TCells = this->T_.internalField();
    const scalarField& pCells = this->p_.internalField();
    
    forAll(TCells, celli)
    {
        //TODO: Cell mixture of only gas species here, or apply mu mixture rule
        muCells[celli] = this->cellMixture(celli).mu( pCells[celli], TCells[celli] );
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        tmu().boundaryField()[patchi] = mu
        (
            this->p_.boundaryField()[patchi],
            this->T_.boundaryField()[patchi], 
            patchi
        );
    }

    return tmu;
}

template<class ThermoType, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::heTwophaseMixtureThermo<ThermoType,MixtureType>::mu
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tmu(new scalarField(T.size()));

    scalarField& mu = tmu();

    forAll(T, facei)
    {
        mu[facei] = this->patchFaceMixture(patchi, facei).mu(p[facei], T[facei]);
    }

    return tmu;
}


template<class ThermoType, class MixtureType> 
Foam::tmp<Foam::volScalarField> 
Foam::heTwophaseMixtureThermo<ThermoType,MixtureType>::dQ_phaseChange() const
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
template<class ThermoType, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heTwophaseMixtureThermo<ThermoType,MixtureType>::nearInterface
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
template<class ThermoType, class MixtureType>
Foam::tmp<Foam::surfaceScalarField>
Foam::heTwophaseMixtureThermo<ThermoType,MixtureType>::surfaceTensionForce()
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
template<class ThermoType, class MixtureType>
void Foam::heTwophaseMixtureThermo<ThermoType,MixtureType>::setPtrs
(
    combustionModels::rhoCombustionModel* combustion
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
template<class ThermoType, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heTwophaseMixtureThermo<ThermoType,MixtureType>::getRefinementField() const
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




template<class ThermoType, class MixtureType>
scalar Foam::heTwophaseMixtureThermo<ThermoType,MixtureType>::solve
(
    volScalarField& rho
)
{
    //Solve for reaction rates
    Info<< "Solving combustion" << endl;
    this->rho_ = alphaLiquid_.rhoAlpha() + alphaVapor_.rhoAlpha();
    this->rho_.correctBoundaryConditions();
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
    
    Info<< "Min,max rho = " << Foam::min(this->rho_).value() << ", " 
        << Foam::max(this->rho_).value() << endl;
        
        
    YSum_ = alphaLiquid_.Yp() + alphaVapor_.Yp();
    
    Info<< "Min,Max Ysum = " << Foam::min(YSum_).value() 
        << ", " << Foam::max(YSum_).value() << endl;
    
    // Calculate diffusion coefficient for each phase
    alphaLiquid_.calculateDs( combustionPtr_->turbulence().muEff() );
    alphaVapor_.calculateDs( combustionPtr_->turbulence().muEff() );
       
    // Solve for subspecie transport within each phase
    scalar FoLiq = alphaLiquid_.solveSubSpecies( this->p_, this->T_, phaseChangeModels_);
    scalar FoVap = alphaVapor_.solveSubSpecies( this->p_, this->T_, phaseChangeModels_);

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



template<class ThermoType, class MixtureType>
Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::heTwophaseMixtureThermo<ThermoType,MixtureType>::TSuSp() const
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
    
    tTSuSp.first()() *= rCp();
    tTSuSp.second()() *= rCp();
    
    return tTSuSp;
}



template<class ThermoType, class MixtureType>
void Foam::heTwophaseMixtureThermo<ThermoType,MixtureType>::solveAlphas
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
    word gradScheme("grad(alphaVaporSmooth)");
    const volVectorField gradAlpha(fvc::grad(alphaLiquid_, gradScheme));
    
    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

    // Face unit interface normal
    surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));

    surfaceScalarField phir(phic*(nHatfv & mesh_.Sf()));
        
    surfaceScalarField phiAlphaL("phiAlphaL",phi_);
    //surfaceScalarField phiAlphaV("phiAlphaV",phi_);
    
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
        divPhaseChange_ += pcmI().Vdot("Vapor") + pcmI().Vdot("Liquid"); //used in the pEqn
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
    
    surfaceScalarField rhoVf(fvc::interpolate(alphaVapor_.rho(this->p_,this->T_)));
    surfaceScalarField rhoLf(fvc::interpolate(alphaLiquid_.rho(this->p_,this->T_)));
    
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


template<class ThermoType, class MixtureType>
void Foam::heTwophaseMixtureThermo<ThermoType,MixtureType>::correct()
{
    if (debug)
    {
        Info<<"entering heTwophaseMixtureThermo<ThermoType,MixtureType>::correct()"<<endl;
    }

    calculate();

    if (debug)
    {
        Info<<"exiting heTwophaseMixtureThermo<ThermoType,MixtureType>::correct()"<<endl;
    }
}


template<class ThermoType, class MixtureType>
bool Foam::heTwophaseMixtureThermo<ThermoType,MixtureType>::read()
{
    if (rhoReactionThermo::read())
    {
        MixtureType::read(*this);
        return true;
    }
    else
    {
        return false;
    }
}



template<class ThermoType, class MixtureType>
tmp<volScalarField> 
Foam::heTwophaseMixtureThermo<ThermoType,MixtureType>::kByCp
(
    const volScalarField& alphat //turbulent k/Cp contribution
) const
{
   // kByCp = (alpha1*k1/Cp1 + alpha2*k2/Cp2) + alphat

    return alphaLiquid_*alphaLiquid_.kappa(this->p_,this->T_)/alphaLiquid_.Cp(this->p_,this->T_)
         + alphaVapor_*alphaVapor_.kappa(this->p_,this->T_)/alphaVapor_.Cp(this->p_,this->T_)
         + alphat;
}


template<class ThermoType, class MixtureType>
tmp<volScalarField> Foam::heTwophaseMixtureThermo<ThermoType,MixtureType>::rCp() const
{
    return 1.0 / this->Cp();
}


template<class ThermoType, class MixtureType>
void Foam::heTwophaseMixtureThermo<ThermoType,MixtureType>::setHE()
{
    scalarField& heCells = this->he_.internalField();
    const scalarField& pCells = this->p_.internalField();
    const scalarField& TCells = this->T_.internalField();

    forAll(heCells, celli)
    {
        heCells[celli] =
            this->cellMixture(celli).HE(pCells[celli], TCells[celli]);
    }

    forAll(this->he_.boundaryField(), patchi)
    {
        this->he_.boundaryField()[patchi] == he
        (
            this->p_.boundaryField()[patchi],
            this->T_.boundaryField()[patchi],
            patchi
        );
    }

    this->heBoundaryCorrection(this->he_);

}


// ************************************************************************* //
