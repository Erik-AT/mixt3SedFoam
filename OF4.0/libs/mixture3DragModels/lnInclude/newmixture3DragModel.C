/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "mixture3DragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mixture3DragModel> Foam::mixture3DragModel::New
(
    const immiscibleIncompressibleThreePhaseMixture& mixture,
	const volVectorField& Ur
)
{
    word modelType(mixture.muModel().viscosityProperties().lookup("dargModel"));

    Info << "Selecting mixture drag model for phase "
        << mixture.phase1Name() 
        << ": "
        << modelType << endl;

    mixtureConstructorTable::iterator cstrIter =
        mixtureConstructorTablePtr_->find(modelType);

    if (cstrIter == mixtureConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown mixture drag model type "
            << modelType << endl << endl
            << "Valid mixture drag model types are : " << endl
            << mixtureConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(mixture, Ur);
}


// ************************************************************************* //
