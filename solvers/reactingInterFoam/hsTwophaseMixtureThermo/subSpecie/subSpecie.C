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

#include "subSpecie.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::subSpecie::subSpecie
(
    const word& name,
    const dictionary& subSpecieDict,
    const fvMesh& mesh,
    volScalarField& specie,
    const gasThermoPhysics& specieData,
    label idx,
    dimensionedScalar phaseSc
)
:
    name_(name),
    subSpecieDict_(subSpecieDict),
    Y_(specie),
    Yp_
    (
        IOobject
        (
            "Yp_" + name,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        specie
    ),
    thermo_(specieData),
    idx_(idx),
    rho0_
    (
        subSpecieDict.lookupOrDefault
        (
            "rho0",
            dimensionedScalar("rho0",dimDensity,0.0)
        )
    ),
    Tc_
    (
        subSpecieDict.lookupOrDefault
        (
            "Tc",
            dimensionedScalar("Tc",dimTemperature,300.0)
        )
    ),
    sigma0_
    (
        subSpecieDict.lookupOrDefault
        (
            "sigma0",
            dimensionedScalar("sigma0",dimMass/dimTime/dimTime,0.0)
        )
    ),
    kappa_
    (
        subSpecieDict.lookupOrDefault
        (
            "kappa",
            dimensionedScalar("kappa",dimPower/dimLength/dimTemperature,0.0)
        )
    ),
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
    D0_
    (
        subSpecieDict.lookupOrDefault
        (
            "D0",
            dimensionedScalar("D0",dimArea/dimTime,0.0)
        )
    ),
    Sc_
    (
        subSpecieDict.lookupOrDefault
        (
            "Sc",
            phaseSc
        )
    )
{

    // Set scalars
    sigmaa_ = (subSpecieDict_.found("sigmaa")) ?
             readScalar(subSpecieDict_.lookup("sigmaa")) : 0;

    if( subSpecieDict_.found("transportModel") )
    {
        // This means this subspecie is a liquid, make sure the other inputs
        // were also provided
        nuModel_ = viscosityModel::New
        (
            "nu" + name_, 
            subSpecieDict_, 
            mesh.lookupObject<volVectorField>("U"), 
            mesh.lookupObject<surfaceScalarField>("phi")
        );
        
        if( rho0_.value() == 0.0   || 
            kappa_.value() == 0.0  || 
            Tc_.value() == 0.0     ||
            sigma0_.value() == 0.0 )
        {
            FatalErrorIn
            (
                "subSpecie::subSpecie"
            )   << "Liquid subspecie " << name
                << "\n    requires entries for rho0, kappa, Tc, and sigma0 "
                << exit(FatalError);
        }
    }
    
    Foam::Info<< "Created and linked subSpecie " << name << Foam::endl;
    Foam::Info<< "  idx = " << idx_ << Foam::endl;
    Foam::Info<< "  rho0 = " << rho0_ << Foam::endl;
    Foam::Info<< "  Tc = " << Tc_ << Foam::endl;
    Foam::Info<< "  sigma0 = " << sigma0_ << Foam::endl;
    Foam::Info<< "  sigmaa = " << sigmaa_ << Foam::endl;
    Foam::Info<< "  kappa = " << kappa_ << Foam::endl;
    
    if( D0_.value() > 0.0 )
    {
        Foam::Info<<"  D0 = " << D0_.value() << Foam::endl;
    }
    else
    {
        Foam::Info<<"  Sc = " << Sc_.value() << Foam::endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::subSpecie::correct()
{
    if( subSpecieDict_.found("transportModel") )
    {
        nuModel_->correct();
    }
}

Foam::autoPtr<Foam::subSpecie> Foam::subSpecie::clone() const
{
    notImplemented("subSpecie::clone() const");
    return autoPtr<subSpecie>(NULL);
}



Foam::tmp<Foam::scalarField> Foam::subSpecie::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCp(new scalarField(T.size()));

    scalarField& cp = tCp();

    forAll(T, facei)
    {
        cp[facei] = thermo_.Cp( T[facei] );
    }

    return tCp;
}

tmp<volScalarField> Foam::subSpecie::Cp
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
                T.mesh().time().timeName(),
                T.mesh()
            ),
            T.mesh(),
            dimensionedScalar("tCp"+name_,dimEnergy/dimMass/dimTemperature,0.0)
        )
    );        
        
    volScalarField& cp = tCp();

    scalarField& cpCells = cp.internalField();
    const scalarField& TCells = T.internalField();

    forAll(TCells, celli)
    {
        cpCells[celli] = thermo_.Cp(TCells[celli]);
    }

    forAll(T.boundaryField(), patchi)
    {
        cp.boundaryField()[patchi] = Cp(T.boundaryField()[patchi], patchi);
    }
        
        
    return tCp;
}


Foam::tmp<Foam::scalarField> Foam::subSpecie::Cv
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCv(new scalarField(T.size()));

    scalarField& cv = tCv();

    forAll(T, facei)
    {
        cv[facei] = thermo_.Cv( T[facei] );
    }

    return tCv;
}

tmp<volScalarField> Foam::subSpecie::Cv
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
                T.mesh().time().timeName(),
                T.mesh()
            ),
            T.mesh(),
            dimensionedScalar("tCv"+name_,dimEnergy/dimMass/dimTemperature,0.0)
        )
    );        
        
    volScalarField& cv = tCv();

    scalarField& cvCells = cv.internalField();
    const scalarField& TCells = T.internalField();

    forAll(TCells, celli)
    {
        cvCells[celli] = thermo_.Cv(TCells[celli]);
    }

    forAll(T.boundaryField(), patchi)
    {
        cv.boundaryField()[patchi] = Cv(T.boundaryField()[patchi], patchi);
    }
        
        
    return tCv;
}

Foam::tmp<Foam::scalarField> Foam::subSpecie::kappa
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tkappa(new scalarField(T.size()));

    forAll(T, facei)
    {
        tkappa()[facei] = thermo_.kappa( T[facei] );
    }

    return tkappa;
}

tmp<volScalarField> Foam::subSpecie::kappa
(
    const volScalarField& T
) const
{
    tmp<volScalarField> tkappa
    (
        new volScalarField
        (
            IOobject
            (
                "tkappa"+name_,
                T.mesh().time().timeName(),
                T.mesh()
            ),
            T.mesh(),
            dimensionedScalar
            (
                "tkappa"+name_,
                dimPower/dimLength/dimTemperature,
                0.0
            )
        )
    );        
        

    scalarField& kappaCells = tkappa().internalField();
    const scalarField& TCells = T.internalField();

    forAll(TCells, celli)
    {
        kappaCells[celli] = thermo_.kappa(TCells[celli]);
    }

    forAll(T.boundaryField(), patchi)
    {
        tkappa().boundaryField()[patchi] = 
            kappa(T.boundaryField()[patchi], patchi);
    }
        
        
    return tkappa;
}

tmp<volScalarField> Foam::subSpecie::sigma
(
    const volScalarField& T
) const
{
    if( sigmaa_ > 0.0 )
    {
        return sigma0_ * Foam::pow(1.0 - Foam::min(T,Tc_)/Tc_,sigmaa_);
    }
    else
    {
        return sigma0_*neg(T-Tc_);
    }
}



tmp<volScalarField> Foam::subSpecie::mu
(
    const volScalarField& T
) const
{
    tmp<volScalarField> tmu
    (
        new volScalarField
        (
            IOobject
            (
                "tmu"+name_,
                T.mesh().time().timeName(),
                T.mesh()
            ),
            T.mesh(),
            dimensionedScalar
            (
                "tmu"+name_,
                dimMass/dimLength/dimTime,
                0.0
            )
        )
    );        
        

    scalarField& muCells = tmu().internalField();
    const scalarField& TCells = T.internalField();

    forAll(TCells, celli)
    {
        muCells[celli] = thermo_.mu(TCells[celli]);
    }

    forAll(T.boundaryField(), patchi)
    {
        tmu().boundaryField()[patchi] = 
            mu(T.boundaryField()[patchi], patchi);
    }
        
        
    return tmu;
}

// Calculate mu on patches
Foam::tmp<Foam::scalarField> Foam::subSpecie::mu
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tmu(new scalarField(T.size()));

    forAll(T, facei)
    {
        tmu()[facei] = thermo_.mu( T[facei] );
    }

    return tmu;
}


void Foam::subSpecie::calculateDs
(
    const volScalarField& mut,
    const volScalarField& rho,
    const volScalarField& T
)
{
    if( D0_.value() > 0.0 )
    {
        D_ = D0_*fvc::interpolate(rho);
    }
    else
    {
        //D = (mu + mut) / Sc
        volScalarField muEff = mut;

        if( hasNuModel() )
        {
            // Liquid is laminar, do not include mut in liquid diffusion
            muEff = nuModel().nu()*rho;
        }
        else
        {
            // Calculate mu using Sutherland transport
            scalarField& muCells = muEff.internalField();
            const scalarField& TCells = T.internalField();

            forAll(TCells, celli)
            {
                muCells[celli] += thermo_.mu(TCells[celli]);
            }

            forAll(T.boundaryField(), patchi)
            {
                muEff.boundaryField()[patchi] += 
                    mu(T.boundaryField()[patchi], patchi);
            }
        }
        
        D_ = fvc::interpolate(muEff) / Sc_;
    }
}


// ************************************************************************* //
