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
    Foam::mixture3DragModels::Wallis

Description
	Drag coefficient based on the drift flux closure. The emperical formula for exponent n 
	is given by Wallis
	Wallis, G.B., One-dimensional two-phase flow. 1969

SourceFiles
    Wallis.C

\*---------------------------------------------------------------------------*/

#ifndef Wallis_H
#define Wallis_H
#include "mixture3DragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mixture3DragModels
{

/*---------------------------------------------------------------------------*\
                           Class Wallis Declaration
\*---------------------------------------------------------------------------*/

class Wallis
:
    public mixture3DragModel
{


	// private data

		//Dispersed phase fraction
		const volScalarField& alphad_;

		//Drag coefficient
		volScalarField CD_;


	//private member functions
		// Calculate and return the drag coefficient
        tmp<volScalarField> calcCD() const;

public:

    //- Runtime type information
    TypeName("Wallis");


    // Constructors

        //- Construct from components
        Wallis
        (
            const immiscibleIncompressibleThreePhaseMixture& mixture,
			const volVectorField& Ur
        );


    //- Destructor
    virtual ~Wallis();


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
