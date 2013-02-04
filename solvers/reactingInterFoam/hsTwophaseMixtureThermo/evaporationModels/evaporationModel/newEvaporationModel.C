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

#include "evaporationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::evaporationModel>
Foam::evaporationModel::New
(
    dictionary specieDict,
    const volScalarField& p,
    const volScalarField& T,
    const phase& alphaL,
    const phase& alphaV
)
{
    word evaporationModelTypeName
    (
        specieDict.lookup("evaporationModel")
    );

    Info<< "Selecting evaporation model "
        << evaporationModelTypeName << endl;

    componentsConstructorTable::iterator cstrIter =
        componentsConstructorTablePtr_
            ->find(evaporationModelTypeName);

    if (cstrIter == componentsConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "evaporationModel::New"
        )   << "Unknown evaporationModel type "
            << evaporationModelTypeName << endl << endl
            << "Valid  evaporationModel are : " << endl
            << componentsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<evaporationModel>(cstrIter()(specieDict, p, T, alphaL, alphaV));
}


// ************************************************************************* //
