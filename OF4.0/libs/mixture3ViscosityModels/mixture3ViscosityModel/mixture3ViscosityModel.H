/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2015 OpenFOAM Foundation
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
    Foam::mixture3ViscosityModels

Description
    A namespace for incompressible mixture3ViscosityModel implementations used with 3-phase model.

Class
    Foam::mixture3ViscosityModel

Description
    An abstract base class for incompressible mixture3ViscosityModels used with 3-phase model.

SourceFiles
    mixture3ViscosityModel.C
    mixture3ViscosityModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef mixture3ViscosityModel_H
#define mixture3ViscosityModel_H

#include "dictionary.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedScalar.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class mixture3ViscosityModel Declaration
\*---------------------------------------------------------------------------*/

class mixture3ViscosityModel
{

protected:

    // Protected data

        word name_;
        dictionary viscosityProperties_;
        
        //Pass U to give access to mesh object which in trun give access to the registerObject
        const volVectorField& U_;

    // Private Member Functions

        //- Disallow copy construct
        mixture3ViscosityModel(const mixture3ViscosityModel&);

        //- Disallow default bitwise assignment
        void operator=(const mixture3ViscosityModel&);


public:

    //- Runtime type information
    TypeName("mixture3ViscosityModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            mixture3ViscosityModel,
            dictionary,
            (
                const word& name,
                const dictionary& viscosityProperties,
                const volVectorField& U
            ),
            (name, viscosityProperties, U)
        );


    // Selectors

        //- Return a reference to the selected viscosity model
        static autoPtr<mixture3ViscosityModel> New
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U
        );


    // Constructors

        //- Construct from components
        mixture3ViscosityModel
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U
        );


    //- Destructor
    virtual ~mixture3ViscosityModel()
    {}


    // Member Functions

        //- Return the phase transport properties dictionary
        const dictionary& viscosityProperties() const
        {
            return viscosityProperties_;
        }


        //- Return the mixture kinematic viscosity given the kinematic viscosity of the continuous  
        //  phase and the packing limit
		
        virtual tmp<volScalarField> nu
		(
			const volScalarField& nuc, 
			const scalar& alphaMax
		) const = 0;

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
