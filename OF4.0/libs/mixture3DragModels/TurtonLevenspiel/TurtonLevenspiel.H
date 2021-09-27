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
    Foam::mixture3DragModels::TurtonLevenspiel

Description
	Turton-Levenspiel drag coefficient
	Turton, R. and Levenspiel, O. (1986). A short note on the drag correlation for spheres. 
	Powder Technology, 47(1), pp.83-86.

SourceFiles
    Turton.C

\*---------------------------------------------------------------------------*/

#ifndef TurtonLevenspiel_H
#define TurtonLevenspiel_H

#include "mixture3DragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mixture3DragModels
{

/*---------------------------------------------------------------------------*\
                           Class TurtonLevenspiel Declaration
\*---------------------------------------------------------------------------*/

class TurtonLevenspiel
:
    public mixture3DragModel
{


	// private data
		//Drag coefficient
		volScalarField CD_;

	//private member functions
		// Calculate and return the drag coefficient
        tmp<volScalarField> calcCD() const;

public:

    //- Runtime type information
    TypeName("TurtonLevenspiel");


    // Constructors

        //- Construct from components
        TurtonLevenspiel
        (
            const immiscibleIncompressibleThreePhaseMixture& mixture,
			const volVectorField& Ur
        );


    //- Destructor
    virtual ~TurtonLevenspiel();


    // Member Functions
		//Return drag coefficient
        tmp<volScalarField> CD() const
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
