/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "threePhaseMixture.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::threePhaseMixture::threePhaseMixture
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    phase1Name_(wordList(dict.lookup("phases"))[0]),
    phase2Name_(wordList(dict.lookup("phases"))[1]),
	phase3Name_(wordList(dict.lookup("phases"))[2]),

    alphad_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase1Name_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    alphac_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase2Name_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    gamma_
    (
        IOobject
        (
            "gamma",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphad_ + alphac_
    ),
    
    

    alphag_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase3Name_),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        1.0 - alphad_ - alphac_,
		alphac_.boundaryField().types()
    )
{}


// ************************************************************************* //
