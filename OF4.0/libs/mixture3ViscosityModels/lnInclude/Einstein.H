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

Class
    Foam::mixture3ViscosityModels::Einstein

Description


SourceFiles
    Einstein.C

\*---------------------------------------------------------------------------*/

#ifndef Einstein_H
#define Einstein_H

#include "mixture3ViscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class incompressibleThreePhaseInteractingMixture;

namespace mixture3ViscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class Einstein Declaration
\*---------------------------------------------------------------------------*/

class Einstein
:
    public mixture3ViscosityModel
{

protected:

    // Protected data

        //- Dispersed phase fraction
        const volScalarField& alpha_;
        

public:

    //- Runtime type information
    TypeName("Einstein");


    // Constructors

        //- Construct from components
        Einstein
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U
        );


    //- Destructor
    ~Einstein()
    {}


    // Member Functions

        //- Return the mixture kinematic viscosity given the kinematic viscosity of the continuous  
        //  phase and the packing limit
		tmp<volScalarField> nu
		(
			const volScalarField& nuc, 
			const scalar& alphaMax
		) const;

        //- Read transportProperties dictionary
        bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mixture3ViscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
