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
    Foam::subSpecie

Description
    Single compressible phase derived from the phase-fraction.
    Used as part of the hsMultiphaseMixtureThermo for interface-capturing
    simulations in conjunction with a reaction subclass. The phase contains
    its constituitive elements (subspecies).

SourceFiles
    subSpecie.C

\*---------------------------------------------------------------------------*/

#ifndef subspecie_H
#define subspecie_H

#include "volFields.H"
#include "PtrDictionary.H"
#include "dictionaryEntry.H"
#include "fvCFD.H"
#include "rhoChemistryCombustionModel.H"
#include "turbulenceModel.H"
#include "incompressible/viscosityModels/viscosityModel/viscosityModel.H"
#include "reactingMixture.H"
#include "reactionTypes.H"
#include "dimensionedScalar.H"
#include "dimensionedScalarFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    //Forward declaration of evaporation model class
    class evaporationModel;
    
/*---------------------------------------------------------------------------*\
                           Class phase Declaration
\*---------------------------------------------------------------------------*/

class subSpecie
{
    // Private data

        //- Subspecie name
        word name_;
        
        //- Initial concentration
        scalar Yinit_;
        
        //- Subspecie dictionary (in thermophysicalProperties)
        dictionary subSpecieDict_;
                       
        //- Viscosity model for this phase
        autoPtr<viscosityModel> nuModel_;
        
        //- Mass fraction of this subspecie
        volScalarField& Y_;
        
        volScalarField Yp_;
        
        //- Thermo properties of this subspecie
        const gasThermoPhysics& thermo_;
        
        //- Specie's position in the global specie list
        label idx_;
        
        dimensionedScalar rho0_;
        dimensionedScalar Tc_;
        dimensionedScalar sigma0_;
        scalar sigmaa_;
        
        scalar sigmaLJ_;
        scalar TLJ_;
        
        dimensionedScalar kappa_;
        
public:

    // Constructors

        //- Construct from components
        subSpecie
        (
            const word& name,
            const dictionary& subSpecieDict,
            const fvMesh& mesh,
            volScalarField& species,
            const gasThermoPhysics& speciesData,
            label idx
        );

        //- Return clone
        autoPtr<subSpecie> clone() const;

        //- Return a pointer to a new phase created on freestore
        //  from Istream
        class iNew
        {
            const fvMesh& mesh_;
            PtrList<volScalarField>& species_;
            const PtrList<gasThermoPhysics>& speciesData_;

        public:

            iNew
            (
                const fvMesh& mesh,
                PtrList<volScalarField>& species,
                const PtrList<gasThermoPhysics>& speciesData
            )
            :
                mesh_(mesh),
                species_(species),
                speciesData_(speciesData)
            {}

            autoPtr<subSpecie> operator()(Istream& is) const
            {
                dictionaryEntry ent(dictionary::null, is);
                label s = -1;
                forAll(species_, i)
                {
                    if( ent.keyword() == species_[i].name() )
                    {
                        s = i;
                        break;
                    }
                }
    
                if( s == -1 )
                {
                    FatalErrorIn
                    (
                        "subSpecie::subSpecie"
                    )   << "Could not match subSpecie " << ent.keyword()
                        << "\n    with any specie in thermo.composition() "
                        << exit(FatalError);
                }
                return autoPtr<subSpecie>
                (
                    new subSpecie
                    (
                        ent.keyword(), 
                        ent, 
                        mesh_, 
                        species_[s], 
                        speciesData_[s], 
                        s
                    )
                );
            }
        };


    // Member Functions

        const dictionary& dict() const
        {
            return subSpecieDict_;
        }
        
        const word& name() const
        {
            return name_;
        }
        
        const word& keyword() const
        {
            return name();
        }
                
        const dimensionedScalar rho0() const
        {
            return rho0_;
        }
        
        const dimensionedScalar sigma0() const
        {
            return sigma0_;
        }
        
        const volScalarField& Y() const
        {
            return Y_;
        }
                
        volScalarField& Y()
        {
            return Y_;
        }
        
        const volScalarField& Yp() const
        {
            return Yp_;
        }
                
        volScalarField& Yp()
        {
            return Yp_;
        }
                
        bool hasNuModel() const
        {
            return subSpecieDict_.found("transportModel");
        }
        
        const viscosityModel& nuModel() const
        {
            return nuModel_();
        }
            
        tmp<volScalarField> sigma(const volScalarField& T) const;
                
        const gasThermoPhysics& thermo() const
        {
            return thermo_;
        }
        
        dimensionedScalar W() const
        {
            return dimensionedScalar("W_"+name_, dimMass/dimMoles, thermo_.W());
        }
        
        const dimensionedScalar& kappa() const
        {
            return kappa_;
        }
        
        label idx() const
        {
            return idx_;
        }
        
        dimensionedScalar RR() const
        {
            return dimensionedScalar("R", dimensionSet(1, 2, -2, -1, -1), thermo_.RR);
        }
        
        tmp<scalarField> Cv
        (
            const scalarField& T,
            const label patchi
        ) const;

        /*tmp<volScalarField> hs
        (
            const volScalarField& T
        ) const;*/
        
        tmp<volScalarField> Cv(const volScalarField& T) const;
        
        void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
