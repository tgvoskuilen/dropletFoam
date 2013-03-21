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
    Foam::hsTwophaseMixtureThermo

Description
    A multiphase reaction thermo class. The various parts and pointers are
    set in the following sequence (after construction, you must call setPtr
    since the thermo class does not have access to the combustion model at
    creation)
    
        hsTwophaseMixtureThermo constructor
            combustionPtr set to NULL (can't access yet)
            phase constructor
                evapPtr set to NULL (needs all phases created first)
                combustionPtr set to NULL (can't access yet)
            for each phase with an evaporation model, create the model and
                set phase.evapPtr (remains NULL for phases without evaporation)
                
        hsTwophaseMixtureThermo.setPtr( combustion )
            sets combustionPtr in thermo and in all its phases
            
            
SourceFiles
    hsTwophaseMixtureThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hsTwophaseMixtureThermo_H
#define hsTwophaseMixtureThermo_H

#include "hsReactionThermo.H"
#include "phase.H"
#include "subSpecie.H"
#include "evaporationModel.H"
#include "PtrDictionary.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "rhoChemistryCombustionModel.H"
#include "multivariateSurfaceInterpolationScheme.H"
#include "subCycle.H"

#include "Time.H"
#include "MULES.H"
#include "fvcSnGrad.H"
#include "fvcFlux.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class hsTwophaseMixtureThermo Declaration
\*---------------------------------------------------------------------------*/

template<class MixtureType>
class hsTwophaseMixtureThermo
:
    public hsReactionThermo,
    public MixtureType
{
    
private:
    // Private Attributes
    
        //- Reference to the host mesh
        const fvMesh& mesh_;
    
        //- Pointer to the combustion model (given after construction)
        combustionModels::rhoChemistryCombustionModel* combustionPtr_;

        //- Vapor phase
        phase alphaVapor_;
        
        //- Liquid phase
        phase alphaLiquid_;
        
        //- Total mass flux
        surfaceScalarField rhoPhi_;
                
        //- Total volume flux
        const surfaceScalarField& phi_;
        
        //- Velocity field
        const volVectorField& U_;
                
        //- Sharpened vapor fraction (for calculating curvature kappaI)
        volScalarField alphaVaporSmooth_;
        
        //- Curvature
        volScalarField kappaI_;
        
        //- Surface tension
        volScalarField sigma_;
                
        //- Evaporation mask to stop evaporation at diffuse liquid interfaces
        volScalarField evap_mask_;
        
        //- Sum of all mass fractions
        volScalarField YSum_;

        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;
        
        scalar phaseRelaxTime_;
        
        scalar phaseMaskTol_;
        
    // Private methods

        //- Calculate T(hs), psi(p,T), mu(p,T), alpha(p,T), rho(p,T), alphas
        void calculate();

        //- Construct as copy (not implemented)
        hsTwophaseMixtureThermo(const hsTwophaseMixtureThermo<MixtureType>&);

        //- Calculate evaporation rates for each evaporation pair
        void calcEvaporation();

        //- Solve the alpha transport equations
        void solveAlphas
        (
            const scalar cAlpha, 
            const label nAlphaCorr, 
            scalar f = 1.0
        );
                
public:

    //- Runtime type information
    TypeName("hsTwophaseMixtureThermo");


    // Constructors

        //- Construct from mesh
        hsTwophaseMixtureThermo(const fvMesh&);


    //- Destructor
    virtual ~hsTwophaseMixtureThermo();


    // Member functions

        //- Return the compostion of the multi-component mixture
        virtual basicMultiComponentMixture& composition()
        {
            return *this;
        }

        //- Return the compostion of the multi-component mixture
        virtual const basicMultiComponentMixture& composition() const
        {
            return *this;
        }

        //- Update properties
        virtual void correct();


        // Fields derived from thermodynamic state variables - These use
        //  a cellMixture based on mass fractions (Ys)

            //- Chemical enthalpy [J/kg]
            virtual tmp<volScalarField> hc() const;


            //- Sensible enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> hs
            (
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Sensible enthalpy for patch [J/kg]
            virtual tmp<scalarField> hs
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cp() const;
                        
            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cv() const;
            
            //- 
            virtual tmp<scalarField> mu
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- 
            virtual tmp<volScalarField> muV() const;
            
            //- 
            virtual tmp<scalarField> kappa
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- 
            virtual tmp<volScalarField> kappaV() const;
            
            
        //- Read thermophysicalProperties dictionary
        virtual bool read();
        
        
        // Multiphase specific methods
            const volScalarField& alphaLiquid() const
            {
                return alphaLiquid_;
            }
            
            const volScalarField& alphaVapor() const
            {
                return alphaVapor_;
            }
            
            volScalarField& T()
            {
                return T_;
            }
                                    
            //- Get mass flux field
            const surfaceScalarField& rhoPhi() const
            {
                return rhoPhi_;
            }

            //- Solve reaction chemistry, evaporation, and phase transport       
            scalar solve(volScalarField& rho);
            
            //- Calculate the surface tension force field
            tmp<surfaceScalarField> surfaceTensionForce();
            
            tmp<volScalarField> W() const;
            
            //- Set the combustion model and evaporation pointers
            void setPtrs
            (
                combustionModels::rhoChemistryCombustionModel* combustion
            );

            //- Get a 0-1 mask for cells near the interface
            tmp<volScalarField> nearInterface
            (
                double lower=0.01, 
                double upper=0.99
            ) const;
            
            //- Get the mesh refinement criteria field
            tmp<volScalarField> getRefinementField() const;
            
            //- Get the total volume source (for pEqn and deltaT)
            tmp<volScalarField> S_evap() const;
            
            Pair<tmp<volScalarField> > SuSp_evap() const;
                                    
            //- Get the total reaction and evaporation enthalpy source
            tmp<volScalarField> Sh() const;
            tmp<volScalarField> dQ_evap() const;
            
            tmp<volScalarField> kByCv(const volScalarField& muEff) const;
            tmp<volScalarField> rCv() const;
            
            void setHs();
            
            //- Calculate the sharpened volume fraction for the vapor phase
            void calculateAlphaVapor();
            
            //- Calculate the surface tension force from the curvature and
            //  surface tension coefficients
            void calculateSurfaceTension();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "hsTwophaseMixtureThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
