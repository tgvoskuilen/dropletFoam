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

#include "LennardJones.H"
#include "addToRunTimeSelectionTable.H"
#include "evaporationModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusionModels
{
    defineTypeNameAndDebug(LennardJones, 0);
    addToRunTimeSelectionTable(diffusionModel, LennardJones, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::diffusionModels::LennardJones::LennardJones
(
    dictionary specieDict
)
:
    diffusionModel(typeName, specieDict),
    sigma_(readScalar(diffusionDict_.lookup("sigma"))),
    T_(readScalar(diffusionDict_.lookup("T")))
{}
    
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField> Foam::diffusionModels::LennardJones::Dij
(
    const subSpecie& sI,
    const subSpecie& sJ
) const
{
    const fvMesh& mesh = sI.Y().mesh();
    const volScalarField& T = mesh.lookupObject<volScalarField>("T");
    const volScalarField& p = mesh.lookupObject<volScalarField>("p");
    
    //TODO: Check that sI and sJ have the same diffusion model
    const LennardJones& dmJ = dynamic_cast<const LennardJones&>(sJ.diffModel());
    
    tmp<volScalarField> Ts = T / 
            dimensionedScalar
            (
                "Tij",
                dimTemperature,
                Foam::sqrt(T_*dmJ.T())
            );
            
    scalar sigmaijSqr = Foam::pow(0.5*(sigma_ + dmJ.sigma()),2.0);
    scalar Wij = Foam::sqrt(0.5*(1.0/sI.W().value() + 1.0/sJ.W().value()));
    
    tmp<volScalarField> OmegaD = Foam::pow(Ts, -0.145) + Foam::pow(0.5+Ts,-2.0);
    
    /*tmp<volScalarField> OmegaD = 
          1.06036 * Foam::pow(Ts, -0.15610)
        + 0.19300 * Foam::exp(-0.47635*Ts)
        + 1.03587 * Foam::exp(-1.52996*Ts)
        + 1.76474 * Foam::exp(-3.89411*Ts);*/ // Neufeld 1972 form
        
        
    dimensionedScalar A("factor",dimensionSet(1, 1, -3, -1.5, 0),0.026693);
        
    return A*Foam::pow(T,1.5)*Wij / (p*sigmaijSqr*OmegaD);
}

// ************************************************************************* //
