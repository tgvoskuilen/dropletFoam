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
    alphaLiquid_.setPhaseMasks(phaseMaskTol_, p_, T_, phaseChangeModels_, area_ );
    alphaVapor_.setPhaseMasks(phaseMaskTol_, p_, T_, phaseChangeModels_, area_ );
    
    // uses mask to set rhoAlpha
    alphaVapor_.correct(p_,T_);
    alphaLiquid_.correct(p_,T_);
    
    // Could use Coutier-Delgosha method for muT from:
    //   O. Coutier-Delgosha, R. Fortes-Patella, and J. L. Reboud, 
    //   “Evaluation of the turbulence model influence on the numerical 
    //   simulations of unsteady cavitation,” 
    //   Journal of Fluids Engineering, vol. 125, no. 1, pp. 38–45, 2003
    //
    mu_ = alphaVapor_.sharp(0.0)*alphaVapor_.mu(p_, T_) +
          alphaLiquid_.sharp(0.0)*alphaLiquid_.mu(p_, T_);
    mu_.correctBoundaryConditions();
    muAll_ = mu_;
    
    psi_ = alphaVapor_.sharp(0.0)*alphaVapor_.cellMask()*alphaVapor_.psi(T_);
    psi_.correctBoundaryConditions();
    
    rho_ = alphaLiquid_.sharp(0.0)*alphaLiquid_.rho(p_,T_) + 
            alphaVapor_.sharp(0.0)*alphaVapor_.rho(p_,T_);
    //rho_ = alphaLiquid_.rhoAlpha() + alphaVapor_.rhoAlpha(); //same as above but includes cellMasks
    rho_.correctBoundaryConditions();
}

template<class MixtureType>
void Foam::hsTwophaseMixtureThermo<MixtureType>::correctInterface()
{
    //Update kappa, sigma_, alphaSmooth, area_
    
    
    //--Generate smooth volume fraction field-----------------------------------
    
    // Step 1 - Sharpen alphaVapor to remove spurious areas and identify the 
    //          interface (consequence, areas with 
    //          alphaVapor < smootherSharpening/2 will not have surface tension)
    alphaVaporSmooth_ = alphaVapor_.sharp(smootherSharpening_);

    // Step 2 - Diffusively smooth the sharpened field
    
    // Target Fo number for diffusive smoothing. This, in conjunction with the
    // number of iterations (nSmootherIters) will diffusive the interface by an 
    // approximate factor of the smallest cell size. It is assumed that the
    // phase interface is approximately 2-3 times the smallest cell size. If you
    // have significantly smaller cells somewhere else in your mesh, you will
    // have to adjust the number of iterations taken.
    scalar Fo = 0.4; //must be less than or equal to 0.5
    
    // The mesh scale stays constant unless the mesh changes       
    dimensionedScalar Deff = Fo / meshArDelta_;
    
    for( label i = 0; i < nSmootherIters_; ++i )
    {
        alphaVaporSmooth_ += Deff * fvc::laplacian(alphaVaporSmooth_);
    }
    alphaVaporSmooth_.min(1.0);
    alphaVaporSmooth_.max(0.0);
    alphaVaporSmooth_.correctBoundaryConditions();
    
    //--Calculate interface curvature field and normal vector-------------------
    
    // Cell gradient of alphaVaporSmooth
    const volVectorField gradAlpha(fvc::grad(alphaVaporSmooth_));

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

    // Face unit interface normal
    surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));
    
    // this is temporary, purely for postprocessing
    interfaceNormal_ = gradAlpha/(mag(gradAlpha) + deltaN_);
    
    // Store interface normal on faces
    nHatf_ = nHatfv & mesh_.Sf();

    // Simple expression for curvature
    kappaI_ = -fvc::div(nHatf_);


    //--Calculate the surface tension field-------------------------------------
    dimensionedScalar one("one",dimLength,1.0);
    tmp<volScalarField> mask = pos(Foam::mag(kappaI_)*one - 1e-4)
                             + neg(alphaVapor_ - 0.1);
    sigma_ = alphaLiquid_.sigma(T_, mask());
        
    
    //--Calculate the phase area------------------------------------------------
    // Hardt method (J. Comp. Phys.)
    tmp<volScalarField> C = alphaLiquid_.sharp(0.0,0.01);
    word gradScheme("grad(alphaVaporSmooth)");
    
    dimensionedScalar eps("eps",dimArea,SMALL);
    tmp<volScalarField> Cp = Foam::mag(fvc::grad(C(),gradScheme));
    
    dimensionedScalar N = fvc::domainIntegrate(Cp()) 
        / (fvc::domainIntegrate((1-C())*(1-C())*Cp()) + eps);

    area_ = N*(1-C())*(1-C())*Cp; 

    // Clip very small areas
    const volScalarField::DimensionedInternalField& V = mesh_.V();
    area_.dimensionedInternalField() *= pos
    (
        area_.dimensionedInternalField()
      - Foam::pow(V,-1.0/3.0)/100.0
    );
    area_.correctBoundaryConditions();
    
    
}

template<class MixtureType>
void Foam::hsTwophaseMixtureThermo<MixtureType>::calcPhaseChange()
{
    Foam::Info << "Calculating phase change zones" << Foam::endl;
    
    if( alphaLiquid_.subSpecies().size() > 1 )
    {
        /*tmp<volScalarField> toverlap
        (
            new volScalarField
            (
                IOobject
                (
                    "toverlap",
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar("toverlap", dimless, 0.0)
            )
        );
        volScalarField& overlap = toverlap();*/
        dimensionedScalar one("one",dimArea,1.0);
        overlap_ *= 0.0;
        
        forAll(noVaporPairs_, pairI)
        {
            const word& specieA = noVaporPairs_[pairI].first();
            const word& specieB = noVaporPairs_[pairI].second();
            
            const volScalarField& YA = alphaLiquid_.subSpecies()[specieA]->Y();
            const volScalarField& YB = alphaLiquid_.subSpecies()[specieB]->Y();
        
            //overlap += pos(YA-SMALL)*pos(YB-SMALL);
            overlap_ += Foam::mag(fvc::grad(YA))*Foam::mag(fvc::grad(YB))*one;
            
        }
        
        //overlap = neg(overlap);
    
        /*forAllIter
        (
            PtrDictionary<subSpecie>, alphaLiquid_.subSpecies(), specieI
        )
        {
            tgradProd() *= Foam::mag( fvc::grad(specieI().Y())*one );
        }*/
        
        //Info<< "Max grad prod = " << Foam::max(tgradProd()).value() << endl;
        
        // only allow evaporation if gradprod < 1e5 OR alphav > 0.5
        //evap_mask_ = pos(neg(tgradProd() - 100000.0)
        //                 + pos(alphaVapor_ - 0.5) - SMALL);
        //evap_mask_ = pos(overlap + pos(alphaVapor_ - 0.5) - SMALL);
        
        // allow evaporation if overlap < 1e3 OR alphaV > 0.7
        evap_mask_ = pos(neg(overlap_ - 1e5)
                         + pos(alphaVapor_ - 0.5) - SMALL);
        //evap_mask_ = pos(alphaVapor_ + 1.0); //no mask
        
    }
    else
    {
        evap_mask_ = pos(alphaVapor_+1.0); //no multi-liquid mask
    }
    
    // Calculate very smoothed alphaV field
    /*volScalarField alphaVSmooth = alphaV_.sharp(0.01);
    dimensionedScalar dA = pow(min(mesh_.V()),2.0/3.0)/20.0;
    
    for( label i = 0; i < 20; ++i )
    {
        alphaVSmooth += dA * fvc::laplacian(alphaVSmooth);
    }
    alphaVSmooth.min(1.0);
    alphaVSmooth.max(0.0);
    alphaVSmooth.correctBoundaryConditions();
    dimensionedScalar sG("SG",dimless/dimLength,SMALL);
    dimensionedScalar Uo("U0",dimVelocity,1.0);
    
    volVectorField UI = fvc::grad(alphaVSmooth) / (mag(fvc::grad(alphaVSmooth)) + sG)*U0;
    surfaceScalarField phiPhase("phiPhase",fvc::interpolate(UI) & mesh_.Sf());*/
    
    
    forAllIter(PtrDictionary<mixturePhaseChangeModel>, phaseChangeModels_, pcmI)
    {
        pcmI().calculate
        (
            evap_mask_,
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
    interfaceNormal_
    (
        IOobject
        (
            "interfaceNormal",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("interfaceNormal", dimless, vector::zero)
    ),
    nHatf_
    (
        IOobject
        (
            "nHatf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("nHatf", dimArea, 0)
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
    muAll_
    (
        IOobject
        (
            "muAll",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("muAll", dimDensity*dimArea/dimTime, 0.0),
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
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("evap_mask", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    overlap_
    (
        IOobject
        (
            "overlap",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("overlap", dimless, 0.0),
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
    DgradY_
    (
        IOobject
        (
            "DgradY",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("DgradY", dimMass/dimTime, 0.0)
    ),
    deltaN_
    (
        "deltaN",
        1e-8/pow(average(mesh.V()), 1.0/3.0)
    ),
    phaseMaskTol_(lookupOrDefault<scalar>("phaseMaskTol",0.01)),
    meshArDelta_("meshArDelta",dimless/dimArea,1e10),
    nSmootherIters_(lookupOrDefault<label>("nSmootherIters",15)),
    smootherSharpening_(lookupOrDefault<scalar>("smootherSharpening",0.1)),
    phaseClipTol_(lookupOrDefault<scalar>("phaseClipTol",1e-6)),
    noVaporPairs_(lookup("noVaporPairs"))
{
    // Check that the noVaporPairs are all valid
    forAll(noVaporPairs_, pairI)
    {
        const word& specieA = noVaporPairs_[pairI].first();
        const word& specieB = noVaporPairs_[pairI].second();
        
        if( !alphaLiquid_.subSpecies().found(specieA) ||
            !alphaLiquid_.subSpecies().found(specieB) )
        {
            FatalErrorIn
            (
                "hsTwophaseMixtureThermo::hsTwophaseMixtureThermo"
            )   << "Species " << specieA << " and " << specieB
                << "\n    are not both liquid subspecies."
                << exit(FatalError);
        }
    }
    
    updateMeshArDelta();
    
    T_.correctBoundaryConditions();
    p_.correctBoundaryConditions();
    
    divPhaseChange_.oldTime();
    divComp_.oldTime();
    alphaLiquid_.setOtherPhase( &alphaVapor_ );
    alphaVapor_.setOtherPhase( &alphaLiquid_ );

    setHs();
    correctInterface();
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
    volScalarField& rho,
    label PIMPLEcorr
)
{
    // Solve for reaction rates
    //  Create a highly sharpened rhoAlpha set to set Y with for reactions
    //
    Info<< "Solving combustion" << endl;
    
    // Make sharp Y fields for reaction calculations
    
    /*tmp<volScalarField> alphaLS = alphaLiquid_.sharp(0.6,0.2);
    tmp<volScalarField> rhoAlphaLS = alphaLiquid_.rho(p_,T_)*alphaLS();
    tmp<volScalarField> rhoAlphaVS = alphaVapor_.rho(p_,T_)*(1-alphaLS());
    
    rho_ = rhoAlphaLS() + rhoAlphaVS();
    rho_.correctBoundaryConditions();
    
    alphaLiquid_.updateGlobalYs( rhoAlphaLS(), rhoAlphaVS() );
    alphaVapor_.updateGlobalYs( rhoAlphaVS(), rhoAlphaLS() );*/
    
    rho_ = alphaLiquid_.rhoAlpha() + alphaVapor_.rhoAlpha();
    rho_.correctBoundaryConditions();
    alphaLiquid_.updateGlobalYs( alphaLiquid_.rhoAlpha(), alphaVapor_.rhoAlpha() );
    alphaVapor_.updateGlobalYs( alphaVapor_.rhoAlpha(), alphaLiquid_.rhoAlpha() );
    
    
    combustionPtr_->correct();

    // Update global mass fractions based on phase-based mass fractions
    /*rho_ = alphaLiquid_.rhoAlpha() + alphaVapor_.rhoAlpha();
    rho_.correctBoundaryConditions();
    alphaLiquid_.updateGlobalYs( alphaLiquid_.rhoAlpha(), alphaVapor_.rhoAlpha() );
    alphaVapor_.updateGlobalYs( alphaVapor_.rhoAlpha(), alphaLiquid_.rhoAlpha() );*/
    
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
    
    /*if( PIMPLEcorr == 0)
    {
        // update interface
        correctInterface();
    }*/
    
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
    alphaLiquid_.calculateDs( combustionPtr_->turbulence().mut(), p_, T_ );
    alphaVapor_.calculateDs( combustionPtr_->turbulence().mut(), p_, T_ );
       
    // Solve for subspecie transport within each phase
    tmp<surfaceScalarField> DgradYLCv = alphaLiquid_.solveSubSpecies( p_, T_, phaseChangeModels_);
    tmp<surfaceScalarField> DgradYVCv = alphaVapor_.solveSubSpecies( p_, T_, phaseChangeModels_);
    
    DgradY_ = (DgradYLCv*fvc::interpolate(alphaLiquid_) + DgradYVCv*fvc::interpolate(alphaVapor_))*fvc::interpolate(rCv());

    // Update global mass fractions based on phase-based mass fractions
    alphaLiquid_.updateGlobalYs( alphaLiquid_.rhoAlpha(), alphaVapor_.rhoAlpha() );
    alphaVapor_.updateGlobalYs( alphaVapor_.rhoAlpha(), alphaLiquid_.rhoAlpha() );
    
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
Foam::hsTwophaseMixtureThermo<MixtureType>::TSuSp
(
    word mode
) const
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
    
    if( mode == "Cv" )
    {
        tTSuSp.first()() *= rCv();
        tTSuSp.second()() *= rCv();
    }
    else
    {
        tTSuSp.first()() *= rCp();
        tTSuSp.second()() *= rCp();
    }
    
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
    
    //Calculate interface curvature field - For computation of the interface
    // normal, we use the smoothed (quasi-mollified) alpha field.
    // interFoam family of solvers does not update nHat here any more, but they
    // used to. Now it's updated after all the alpha sub-cycles
    //correctInterface();
           
    surfaceScalarField phir(-phic*nHatf_);
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
            divU - divPhaseChange_
        );
        
        volScalarField::DimensionedInternalField Su
        (
            IOobject
            (
                "Su",
                mesh_.time().timeName(),
                mesh_
            ),
            Sv
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
        
    
    surfaceScalarField rhoVf(fvc::interpolate(alphaVapor_.rho(p_,T_)));
    surfaceScalarField rhoLf(fvc::interpolate(alphaLiquid_.rho(p_,T_)));
    
    alphaLiquid_.rhoPhiAlpha() += f * phiAlphaL * rhoLf;
    alphaVapor_.rhoPhiAlpha() += f * (phi_ - phiAlphaL) * rhoVf;

    rhoPhi_ += f * (phiAlphaL*(rhoLf - rhoVf) + phi_*rhoVf);

        
    // Re-sharpen alpha field
    //  This is not mass conserving, but prevents excessive floatsom. If
    //  phaseClipTol is set to 0, this merely keeps alphaL between 0 and 1
    volScalarField& alphaL = alphaLiquid_;
    alphaL = alphaLiquid_.sharp(phaseClipTol_);
    alphaVapor_ == scalar(1) - alphaLiquid_;
        
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
        Info<<"Reading smoother iterations and clipping tolerance" << endl;
        nSmootherIters_     = lookupOrDefault<label>("nSmootherIters",15);
        smootherSharpening_ = lookupOrDefault<scalar>("smootherSharpening",0.1);
        phaseClipTol_       = lookupOrDefault<scalar>("phaseClipTol",1e-6);
        
        MixtureType::read(*this);
        return true;
    }
    else
    {
        return false;
    }
}



template<class MixtureType>
tmp<volScalarField> 
Foam::hsTwophaseMixtureThermo<MixtureType>::kByCp
(
    const volScalarField& alphat //turbulent k/Cp contribution
) const
{
   // kByCp = (alpha1*k1/Cp1 + alpha2*k2/Cp2) + alphat

    //return alphaLiquid_*alphaLiquid_.kappa(T_)/alphaLiquid_.Cp(T_)
    //     + alphaVapor_*alphaVapor_.kappa(T_)/alphaVapor_.Cp(T_)
    //     + alphat;
         
    return rCp()*
           (
               alphaLiquid_*alphaLiquid_.kappa(T_)
             + alphaVapor_*alphaVapor_.kappa(T_)
           )
           + alphat; 
}


template<class MixtureType>
tmp<volScalarField> Foam::hsTwophaseMixtureThermo<MixtureType>::rCp() const
{
    //return 1.0 / Cp();
    //return alphaLiquid_/alphaLiquid_.Cp(T_) + alphaVapor_/alphaVapor_.Cp(T_);
    
    return (alphaLiquid_.rhoAlpha() + alphaVapor_.rhoAlpha()) / 
           (
               alphaLiquid_.rhoAlpha()*alphaLiquid_.Cp(T_)
             + alphaVapor_.rhoAlpha()*alphaVapor_.Cp(T_)
           );
}

template<class MixtureType>
tmp<volScalarField> 
Foam::hsTwophaseMixtureThermo<MixtureType>::kByCv
(
    const volScalarField& alphat //turbulent k/Cv contribution
) const
{
    //return alphaLiquid_*alphaLiquid_.kappa(T_)/alphaLiquid_.Cv(T_)
    //     + alphaVapor_*alphaVapor_.kappa(T_)/alphaVapor_.Cv(T_)
    //     + alphat;
         
    return rCv()*
           (
               alphaLiquid_*alphaLiquid_.kappa(T_)
             + alphaVapor_*alphaVapor_.kappa(T_)
           )
           + alphat; 
}


template<class MixtureType>
tmp<volScalarField> Foam::hsTwophaseMixtureThermo<MixtureType>::rCv() const
{
    //return alphaLiquid_/alphaLiquid_.Cv(T_) + alphaVapor_/alphaVapor_.Cv(T_);
    
    return (alphaLiquid_.rhoAlpha() + alphaVapor_.rhoAlpha()) / 
           (
               alphaLiquid_.rhoAlpha()*alphaLiquid_.Cv(T_)
             + alphaVapor_.rhoAlpha()*alphaVapor_.Cv(T_)
           );
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

template<class MixtureType>
void Foam::hsTwophaseMixtureThermo<MixtureType>::updateMeshArDelta()
{
    // Update the laplacian stability factor for explicit diffusive smoothing
    tmp<fvScalarMatrix> lap = fvm::laplacian(alphaVaporSmooth_);
    
    tmp<volScalarField> sumAnb = lap().H1();
    
    scalar localMax = Foam::max(Foam::mag(sumAnb())).value();
    
    Foam::reduce(localMax, maxOp<scalar>());
    
    meshArDelta_.value() = localMax;
}


// ************************************************************************* //
