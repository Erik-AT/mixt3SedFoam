/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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

#include "KriegerDougherty.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mixture3ViscosityModels
{
    defineTypeNameAndDebug(KriegerDougherty, 0);

    addToRunTimeSelectionTable
    (
        mixture3ViscosityModel,
        KriegerDougherty,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixture3ViscosityModels::KriegerDougherty::KriegerDougherty
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U
)
:
    mixture3ViscosityModel(name, viscosityProperties, U),
    alpha_
    (
        U.mesh().lookupObject<volScalarField>
        (
            IOobject::groupName
            (
                viscosityProperties.lookupOrDefault<word>("alpha", "alpha"),
                viscosityProperties.dictName()
            )
        )
    )    
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::mixture3ViscosityModels::KriegerDougherty::nu
(
	const volScalarField& nuc, 
	const scalar& alphaMax
) const
{
    return
    (
		nuc*pow(1.0 - min(alpha_/alphaMax,scalar(0.999)),-2.5*alphaMax)
    );
}


bool Foam::mixture3ViscosityModels::KriegerDougherty::read
(
    const dictionary& viscosityProperties
)
{
    return true;
}


// ************************************************************************* //
