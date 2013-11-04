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
    Foam::viscosityModels::NewtonianPolynomial

Description
    An incompressible Newtonian viscosity model with exponential dependence 
    on temperature.
    
    Ts = (T-Tb)/(Tc-Tb)
    nu = nuB*(exp(-3.2*Ts)+0.15)
        

SourceFiles
    NewtonianPolynomial.C

\*---------------------------------------------------------------------------*/

#ifndef NewtonianPolynomial_H
#define NewtonianPolynomial_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class Newtonian Declaration
\*---------------------------------------------------------------------------*/

class NewtonianPolynomial
:
    public viscosityModel
{
    // Private data

        dimensionedScalar nuB_;
        dimensionedScalar Tb_;
        dimensionedScalar Tc_;
        
        //nu_ = nuB*(exp(-3.2*Ts)+0.15)

        volScalarField nu_;


public:

    //- Runtime type information
    TypeName("NewtonianPolynomial");


    // Constructors

        //- construct from components
        NewtonianPolynomial
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    ~NewtonianPolynomial()
    {}


    // Member Functions

        //- Return the laminar viscosity
        tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Correct the laminar viscosity (not appropriate, viscosity constant)
        void correct();

        //- Read transportProperties dictionary
        bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
