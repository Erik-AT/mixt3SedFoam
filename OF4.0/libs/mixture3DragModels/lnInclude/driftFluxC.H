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
    Foam::mixture3DragModels::driftFluxC

Description


SourceFiles
    driftFluxC.C

\*---------------------------------------------------------------------------*/

#ifndef driftFluxC_H
#define driftFluxC_H
#include "mixture3DragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mixture3DragModels
{

/*---------------------------------------------------------------------------*\
                           Class driftFluxC Declaration
\*---------------------------------------------------------------------------*/

class driftFluxC
:
    public mixture3DragModel
{


	// private data

		//Dispersed phase fraction
		const volScalarField& alphad_;

		//Model exponent: Default = 2.65, Calibrated for the 1DSedimentation case = 3.15
		scalar m_;

		//Minimum volume fraction to activate frictional pressure
		dimensionedScalar alphaMinFr_;

		scalar alphaMax_;

		//Drag coefficient
		volScalarField CD_;

	//private member functions
		// Calculate and return the drag coefficient
        tmp<volScalarField> calcCD() const;

public:

    //- Runtime type information
    TypeName("driftFluxC");


    // Constructors

        //- Construct from components
        driftFluxC
        (
            const immiscibleIncompressibleThreePhaseMixture& mixture,
			const volVectorField& Ur
        );


    //- Destructor
    virtual ~driftFluxC();


    // Member Functions
		//Return drag coefficient
        virtual tmp<volScalarField> CD() const
		{
			return CD_;
		}

		//Correct (recalculate) drag coefficient
		virtual void correct() 
		{
			CD_ = calcCD();
		}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mixture3DragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
