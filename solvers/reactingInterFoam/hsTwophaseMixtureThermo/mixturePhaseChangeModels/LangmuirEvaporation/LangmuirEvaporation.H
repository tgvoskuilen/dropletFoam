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
    Foam::LangmuirEvaporation

Description
    Hertz Knudsen pressure-driven evaporation model

SourceFiles
    LangmuirEvaporation.C

\*---------------------------------------------------------------------------*/

#ifndef LangmuirEvaporation_H
#define LangmuirEvaporation_H

#include "mixturePhaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mixturePhaseChangeModels
{
/*---------------------------------------------------------------------------*\
                      Class LangmuirEvaporation Declaration
\*---------------------------------------------------------------------------*/

class LangmuirEvaporation
:
    public mixturePhaseChangeModel
{

    //model-specific stuff
        volScalarField x_;
        volScalarField xL_;
        volScalarField coeffC_;
        volScalarField coeffV_;
        volScalarField p_vap_;
                
        //- Vapor pressure coefficients for boiling to critical temperature
        //  range
        dimensionedScalar Pc_;
        dimensionedScalar Tc_;
        dimensionedScalar Tb_;
        dimensionedScalar Lb_;
        scalar La_;
        vector PvCoeffs_;
        scalar betaV_;
        scalar betaC_;
        
        word vapor_specie_;
        word liquid_specie_;
        
        dimensionedScalar W_;

public:

    //- Runtime type information
    TypeName("LangmuirEvaporation");

    // Constructors

        //- construct from components
        LangmuirEvaporation
        (
            const word name,
            const fvMesh& mesh,
            const phase& alphaL,
            const phase& alphaV,
            const PtrList<gasThermoPhysics>& speciesData,
            dictionary phaseChangeDict
        );
        
    //- Destructor
        virtual ~LangmuirEvaporation() {}


    // Member Functions
        virtual void calculate
        (
            const volScalarField& evapMask,
            const volScalarField& area
        );
        
        tmp<volScalarField> dPvdT() const;
        tmp<volScalarField> L() const;
        tmp<volScalarField> dLdT() const;
        
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
