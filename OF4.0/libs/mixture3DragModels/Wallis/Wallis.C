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

#include "Wallis.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mixture3DragModels
{
    defineTypeNameAndDebug(Wallis, 0);

    addToRunTimeSelectionTable
    (
        mixture3DragModel,
        Wallis,
        mixture
    );
}
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> 
Foam::mixture3DragModels::Wallis::calcCD() const
{
	//Single partilce Cd
	Foam::volScalarField Cd0											
							(
								  24.0*(1.0+0.173*pow(Re(), 0.657))
								/ Re()+0.413/(1.0+16300*pow(Re(), -1.09))
							); 

	Foam::volScalarField n
							(
								  4.7
								* (1.0+0.15*pow(Re(), 0.687))
								/ (1.0+0.253*pow(Re(), 0.687))
							);

	return  (Cd0*pow(max(1.0-alphad_,scalar(1.0e-3)), 3.0-2*n));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixture3DragModels::Wallis::Wallis
(
    const immiscibleIncompressibleThreePhaseMixture& mixture,
	const volVectorField& Ur
)
:
    mixture3DragModel(mixture, Ur),
	alphad_(mixture.alphad()),
CD_
(
	IOobject
    (
        "CD",
        Ur_.time().timeName(),
        Ur_.db(),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    calcCD()
)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mixture3DragModels::Wallis::~Wallis()
{}



// ************************************************************************* //
