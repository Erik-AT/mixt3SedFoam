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

#include "driftFluxC.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mixture3DragModels
{
    defineTypeNameAndDebug(driftFluxC, 0);

    addToRunTimeSelectionTable
    (
        mixture3DragModel,
        driftFluxC,
        mixture
    );
}
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> 
Foam::mixture3DragModels::driftFluxC::calcCD() const
{
	//Single partilce Cd
	Foam::volScalarField Cd0											
							(neg(Re() - 1000)*(24.0*(1.0 + 0.15*pow(Re(), 0.687))/Re())
					  		 + pos(Re() - 1000)*0.44
							);

	Foam::volScalarField mExp( m_ * (
										  neg(alphad_ - alphaMinFr_) 
										+ pos(alphad_ - alphaMinFr_) 
										* Foam::exp(
														(alphad_ - alphaMinFr_)
													  / max(alphaMax_ - alphad_, scalar(1.0e-20)) 
												   )  
									)
							  );

	return  (Cd0*pow(scalar(1.0)-alphad_, -mExp)); 

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixture3DragModels::driftFluxC::driftFluxC
(
    const immiscibleIncompressibleThreePhaseMixture& mixture,
	const volVectorField& Ur
)
:
    mixture3DragModel(mixture, Ur),
	alphad_(mixture.alphad()),
m_(readScalar(mixture.muModel().viscosityProperties().lookup("m"))),
alphaMinFr_(mixture.subDict("granularRheology").lookup("alphaMinFriction")),
alphaMax_(readScalar(mixture.muModel().viscosityProperties().lookup("alphaMax"))),
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

Foam::mixture3DragModels::driftFluxC::~driftFluxC()
{}



// ************************************************************************* //
