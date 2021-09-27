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

#include "Andersson.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mixture3DragModels
{
    defineTypeNameAndDebug(Andersson, 0);

    addToRunTimeSelectionTable
    (
        mixture3DragModel,
        Andersson,
        mixture
    );
}
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> 
Foam::mixture3DragModels::Andersson::calcCD() const
{
	//Foam::volScalarField Re(max(mag(Ur_)*dd_/nuc_, scalar(1.0e-3))); 

	Foam::volScalarField Cd0(24.0*(1.0+0.173*pow(Re(), 0.657))/Re()+0.413/(1.0+16300*pow(Re(), -1.09))); //Single partilce Cd

	Foam::volScalarField Psi(0.5+Foam::atan(262.5*(alphad_-scalar(0.05)))/constant::mathematical::pi); //Smoothing function

	Foam::volScalarField CdI((Cd0-24/Re())/8); //inertial part of Cd

	Foam::volScalarField zq2(((1-alphad_)/max(2*alphad_, scalar(1.0e-3)))*Foam::exp(2.5*alphad_/(1-0.61*alphad_)));

	Foam::volScalarField q3(5*pow(max(alphad_, scalar(0))/max((1-alphad_),scalar(1.0e-3)), 0.45));

	Foam::volScalarField Cd(24*zq2*alphad_/(Re()*max(1-alphad_, scalar(1.0e-3)))+8*CdI*q3);

	return  (Psi*Cd+(1-Psi)*Cd0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixture3DragModels::Andersson::Andersson
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

Foam::mixture3DragModels::Andersson::~Andersson()
{}



// ************************************************************************* //
