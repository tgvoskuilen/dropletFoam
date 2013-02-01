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

#include "droplet.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::droplet::droplet
(
    const word& name,
    const dictionary& dropletDict,
    const fvMesh& mesh
)
:
    mesh_(mesh),
    name_(name),
    dropletDict_(dropletDict),
    center_(dropletDict.lookup("center")),
    radius_(dropletDict.lookup("radius")),
    dV_(readScalar(dropletDict.lookup("delVapor"))),
    Uinit_(dropletDict.lookup("U")),
    liquidSpecies_(dropletDict.lookup("liquidSpecies")),
    vaporSpecies_(dropletDict.lookup("vaporSpecies"))
{
    Foam::Info<< "Created droplet " << name << Foam::endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::droplet> Foam::droplet::clone() const
{
    notImplemented("droplet::clone() const");
    return autoPtr<droplet>(NULL);
}

Foam::List<Foam::word> Foam::droplet::species() const
{
    List<word> species = liquidSpecies_.toc();
    species.append(vaporSpecies_.toc());
    return species;
}


Foam::tmp<Foam::scalarField> Foam::droplet::r() const
{
    vector sr = radius_;
        
    const volVectorField& cellCenters = mesh_.C();
    
    //vector Isr(1/(sr.x()+SMALL), 1/(sr.y()+SMALL), 1/(sr.z()+SMALL));
    
    if (sr.x() < SMALL) { sr.x() = GREAT; }
    if (sr.y() < SMALL) { sr.y() = GREAT; }
    if (sr.z() < SMALL) { sr.z() = GREAT; }
    
    tensor Isr(1/sr.x(), 0, 0, 
               0, 1/sr.y(), 0,
               0, 0, 1/sr.z());
    
    tmp<scalarField> r = Foam::mag(Isr & (cellCenters.internalField() - center_));
    
    
    
    return r;
}

Foam::tmp<Foam::scalarField> Foam::droplet::dropMask() const
{    
    return neg(r() - 1.0);
}

Foam::tmp<Foam::scalarField> Foam::droplet::vaporMask() const
{
    scalar d_layer = Foam::mag(radius_)/dV_;

    tmp<scalarField> f = Foam::max(1.0 - d_layer*(r() - 1.0), 0.0);
    
    return f * (1.0 - dropMask());
}


void Foam::droplet::set
(
    Foam::volScalarField& alphaLiquid,
    Foam::volVectorField& U,
    PtrList<volScalarField>& species
) const
{
    alphaLiquid.internalField() += dropMask();
    U.internalField() += Uinit_*dropMask();

    forAll(species, i)
    {
        if( liquidSpecies_.found(species[i].name()) )
        {
            species[i].internalField() += dropMask()*liquidSpecies_[species[i].name()];
        }
        
        if( vaporSpecies_.found(species[i].name()) )
        {
            species[i].internalField() += vaporMask()*vaporSpecies_[species[i].name()];
        }
    }
}



// ************************************************************************* //
