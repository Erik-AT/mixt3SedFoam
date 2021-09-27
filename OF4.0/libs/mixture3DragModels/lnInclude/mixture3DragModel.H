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

Namespace
    Foam::mixture3DragModels

Description
    A namespace for mixture3DragModel implementations used with 3-phase model.

Class
    Foam::mixture3DragModel

Description
    An abstract base class for mixture3DragModels used with 3-phase model.

SourceFiles
    mixture3DragModel.C
    mixture3DragModelNew.C


\*---------------------------------------------------------------------------*/

#ifndef mixture3DragModel_H
#define mixture3DragModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "immiscibleIncompressibleThreePhaseMixture.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class mixture3DragModel Declaration
\*---------------------------------------------------------------------------*/

class mixture3DragModel
{
protected:

    // Protected data

		//dispersed phase diameter
		dimensionedScalar dd_;

		//continuous phase viscosity
		const volScalarField& nuc_;

        //Relative Velocity between continuous phase and dispersed phase
		const volVectorField& Ur_;
 
      //dimensionedScalar residualPhaseFraction_;
      //dimensionedScalar residualSlip_;

public:

    //- Runtime type information
    TypeName("mixture3DragModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            mixture3DragModel,
            mixture,
            (
            const immiscibleIncompressibleThreePhaseMixture& mixture,
			const volVectorField& Ur 
            ),
            (mixture, Ur)
        );


    // Constructors

        mixture3DragModel
        (
            const immiscibleIncompressibleThreePhaseMixture& mixture,
			const volVectorField& Ur 
        );


    //- Destructor
    virtual ~mixture3DragModel();


    // Selectors

        static autoPtr<mixture3DragModel> New
        (
            const immiscibleIncompressibleThreePhaseMixture& mixture,
			const volVectorField& Ur 
        );


    // Member Functions

		//Return drag coefficient
        virtual tmp<volScalarField> CD() const = 0;

		//Correct drag coefficient
		virtual void correct() = 0;

		//Return Reynolds number
		volScalarField Re() const { return (max(mag(Ur_)*dd_/nuc_, scalar(1.0e-3)));}
       
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
