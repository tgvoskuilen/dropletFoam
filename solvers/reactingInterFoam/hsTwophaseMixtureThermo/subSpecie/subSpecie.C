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
    label idx
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
    )
{

    // Set scalars
    sigmaa_ = (subSpecieDict_.found("sigmaa")) ?
             readScalar(subSpecieDict_.lookup("sigmaa")) : 0;

    if( subSpecieDict_.found("transportModel") )
    {
        nuModel_ = viscosityModel::New
        (
            "nu" + name_, 
            subSpecieDict_, 
            mesh.lookupObject<volVectorField>("U"), 
            mesh.lookupObject<surfaceScalarField>("phi")
        );
    }
        
    Foam::Info<< "Created and linked subSpecie " << name << Foam::endl;
    Foam::Info<< "  idx = " << idx_ << Foam::endl;
    Foam::Info<< "  rho0 = " << rho0_ << Foam::endl;
    Foam::Info<< "  Tc = " << Tc_ << Foam::endl;
    Foam::Info<< "  sigma0 = " << sigma0_ << Foam::endl;
    Foam::Info<< "  sigmaa = " << sigmaa_ << Foam::endl;
    Foam::Info<< "  kappa = " << kappa_ << Foam::endl;
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

/*
tmp<volScalarField> Foam::subSpecie::hs
(
    const volScalarField& T
) const
{
    tmp<volScalarField> ths
    (
        new volScalarField
        (
            IOobject
            (
                "ths"+name_,
                T.mesh().time().timeName(),
                T.mesh()
            ),
            T.mesh(),
            dimensionedScalar("ths"+name_, dimEnergy/dimMass, 0.0)
        )
    );

    forAll(ths(), cellI)
    {
        ths()[cellI] = thermo_.Hs( T[cellI] );
    }
        
    return ths;
}*/

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

Foam::tmp<Foam::scalarField> Foam::subSpecie::kappa
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tkappa(new scalarField(T.size()));

    //scalarField& kappa = tkappa();

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
    return sigma0_ * Foam::pow(1.0 - Foam::min(T,Tc_)/Tc_,sigmaa_);
}

// ************************************************************************* //
