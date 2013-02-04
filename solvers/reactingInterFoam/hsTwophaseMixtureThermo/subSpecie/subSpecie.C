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
#include "evaporationModel.H"
#include "diffusionModel.H"

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
    D_
    (
        IOobject
        (
            "D"+name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("D"+name, dimArea/dimTime, 0.0)
    ),
    thermo_(specieData),
    idx_(idx),
    rho0_(subSpecieDict.lookupOrDefault("rho0",dimensionedScalar("rho0",dimDensity,0.0))),
    Tc_(subSpecieDict.lookupOrDefault("Tc",dimensionedScalar("Tc",dimTemperature,300.0))),
    sigma0_(subSpecieDict.lookupOrDefault("sigma0",dimensionedScalar("sigma0",dimMass/dimTime/dimTime,0.0))),
    kappa_(subSpecieDict.lookupOrDefault("kappa",dimensionedScalar("kappa",dimPower/dimLength/dimTemperature,0.0)))
{

    // Set scalars
    sigmaa_ = (subSpecieDict_.found("sigmaa")) ? readScalar(subSpecieDict_.lookup("sigmaa")) : 0;

    diffusionModel_ = diffusionModel::New(subSpecieDict_);
    
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

tmp<volScalarField> Foam::subSpecie::Dij
(
    const subSpecie& specieJ,
    const compressible::turbulenceModel& turb
) const
{
    return diffusionModel_().Dij( *this, specieJ, turb );
}

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

tmp<volScalarField> Foam::subSpecie::sigma
(
    const volScalarField& T
) const
{
    return sigma0_ * Foam::pow(1.0 - Foam::min(T,Tc_)/Tc_,sigmaa_);
}

// ************************************************************************* //
