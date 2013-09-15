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
    Foam::ThermalDecompReaction

Description
    Thermal decomposition reaction

SourceFiles
    ThermalDecompReaction.C

\*---------------------------------------------------------------------------*/

#ifndef ThermalDecompReaction_H
#define ThermalDecompReaction_H

#include "mixturePhaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mixturePhaseChangeModels
{
/*---------------------------------------------------------------------------*\
                      Class ThermalDecompReaction Declaration
\*---------------------------------------------------------------------------*/

class ThermalDecompReaction
:
    public mixturePhaseChangeModel
{

    //model-specific stuff
        volScalarField xL_;
        volScalarField p_vap_;
                
        //- Vapor pressure coefficients for boiling to critical temperature
        //  range
        dimensionedScalar Pc_;
        dimensionedScalar Tc_;
        dimensionedScalar Tb_;
        dimensionedScalar Lb_;
        dimensionedScalar dH_;
        vector PvCoeffs_;
        dimensionedScalar betaTD_;
        
        //word vapor_specie_;
        word liquid_specie_;
        
        dimensionedScalar W_;

public:

    //- Runtime type information
    TypeName("ThermalDecompReaction");

    // Constructors

        //- construct from components
        ThermalDecompReaction
        (
            const word name,
            const fvMesh& mesh,
            const phase& alphaL,
            const phase& alphaV,
            const PtrList<gasThermoPhysics>& speciesData,
            dictionary phaseChangeDict
        );
        
    //- Destructor
        virtual ~ThermalDecompReaction() {}


    // Member Functions
        virtual void calculate
        (
            const volScalarField& evapMask,
            const volScalarField& area
        );
        
        virtual tmp<volScalarField> mdot(const word& phaseName) const;
        virtual tmp<volScalarField> Vdot(const word& phaseName) const;
        
        virtual Pair<tmp<volScalarField> > YSuSp(const word& specie) const;
        virtual Pair<tmp<volScalarField> > TSuSp() const;
};  


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace mixturePhaseChangeModel

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
