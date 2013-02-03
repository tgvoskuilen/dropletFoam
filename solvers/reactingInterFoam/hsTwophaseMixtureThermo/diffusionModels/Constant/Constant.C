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

#include "Constant.H"
#include "addToRunTimeSelectionTable.H"
#include "evaporationModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusionModels
{
    defineTypeNameAndDebug(Constant, 0);
    addToRunTimeSelectionTable(diffusionModel, Constant, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::diffusionModels::Constant::Constant
(
    dictionary specieDict
)
:
    diffusionModel(typeName, specieDict),
    constD_(diffusionDict_.lookup("D"))
{}
    
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField> Foam::diffusionModels::Constant::Dij
(
    const subSpecie& sI,
    const subSpecie& sJ
) const
{
    const fvMesh& mesh = sI.Y().mesh();
    
    tmp<volScalarField> tDij
    (
        new volScalarField
        (
            IOobject
            (
                "tDij",
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            constD_
        )
    );
    
    return tDij;
}

// ************************************************************************* //
