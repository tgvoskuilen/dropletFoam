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

Class
    Foam::phase

Description
    Single compressible phase derived from the phase-fraction.
    Used as part of the hsMultiphaseMixtureThermo for interface-capturing
    simulations in conjunction with a reaction subclass. The phase contains
    its constituitive elements (subspecies).

SourceFiles
    droplet.C

\*---------------------------------------------------------------------------*/

#ifndef droplet_H
#define droplet_H

#include "volFields.H"
#include "PtrDictionary.H"
#include "dictionaryEntry.H"
#include "fvCFD.H"
#include "dimensionedScalar.H"
#include "dimensionedScalarFwd.H"
#include "cuttableCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class droplet Declaration
\*---------------------------------------------------------------------------*/

class droplet
{
    // Private data
    
        const fvMesh& mesh_;

        //- Droplet name
        word name_;
        
        //- Droplet dictionary
        dictionary dropletDict_;

        vector center_;
        vector radius_;
        scalar dV_;
        vector Uinit_;
        HashTable<scalar> liquidSpecies_;
        HashTable<scalar> vaporSpecies_;
        
        volScalarField dropMask_; //TODO: Change to scalarField
        volScalarField vaporMask_; //TODO: Change to scalarField
        
public:

    // Constructors

        //- Construct from components
        droplet
        (
            const word& name,
            const dictionary& dropletDict,
            const fvMesh& mesh
        );

        //- Return clone
        autoPtr<droplet> clone() const;

        //- Return a pointer to a new droplet created on freestore
        //  from Istream
        class iNew
        {
            const fvMesh& mesh_;

        public:

            iNew
            (
                const fvMesh& mesh
            )
            :
                mesh_(mesh)
            {}

            autoPtr<droplet> operator()(Istream& is) const
            {
                dictionaryEntry ent(dictionary::null, is);
                return autoPtr<droplet>(new droplet(ent.keyword(), ent, mesh_));
            }
        };


    // Member Functions

        const dictionary& dict() const
        {
            return dropletDict_;
        }

        const word& name() const
        {
            return name_;
        }

        const word& keyword() const
        {
            return name();
        }
        
        const volScalarField& dropMask() const
        {
            return dropMask_;
        }
        
        const volScalarField& vaporMask() const
        {
            return vaporMask_;
        }
        
        
        void set
        (
            volScalarField& alphaLiquid,
            volVectorField& U,
            PtrList<volScalarField>& species
        );
        
        void calcDropMask();
        
        List<word> species() const;
        
        tmp<scalarField> r() const;

        void calculate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
