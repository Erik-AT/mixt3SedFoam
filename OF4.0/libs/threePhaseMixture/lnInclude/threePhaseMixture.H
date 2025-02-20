/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::threePhaseMixture

Description
    A three-phase mixture model

SourceFiles
    threePhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef threePhaseMixture_H
#define threePhaseMixture_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class threePhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class threePhaseMixture
{
protected:

    // Protected data

		//Dispersed phase name
        word phase1Name_;

		//Carrier phase name
        word phase2Name_;

		//Gaseous phase name
		word phase3Name_;

		//Dispersed phase volume fraction
        volScalarField alphad_;

        //Carrier phase volume fraction
        volScalarField alphac_;

		//Suspension (carrier + dispersed) volume fraction
        volScalarField gamma_;
        
		//Gaseous phase volume fraction
        volScalarField alphag_;


public:

    // Constructors

        //- Construct from components
        threePhaseMixture
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    ~threePhaseMixture()
    {}


    // Member Functions

        const word& phase1Name() const
        {
            return phase1Name_;
        }

        const word& phase2Name() const
        {
            return phase2Name_;
        }

        const word& phase3Name() const
        {
            return phase3Name_;
        }

        //- Return the phase-fraction of dispersed phase
        const volScalarField& alphad() const
        {
            return alphad_;
        }

        //- Return the phase-fraction of dispersed phase
        volScalarField& alphad()
        {
            return alphad_;
        }

        //- Return the phase-fraction of suspension
        const volScalarField& gamma() const
        {
            return gamma_;
        }

        //- Return the phase-fraction of suspension
        volScalarField& gamma()
        {
            return gamma_;
        }

        //- Return the phase-fraction of carrier phase
        volScalarField& alphac()
        {
            return alphac_;
        }

        //- Return the phase-fraction of carrier phase
        const volScalarField& alphac() const
        {
            return alphac_;
        }

        //- Return the phase-fraction of gaseous phase
        const volScalarField& alphag() const
        {
            return alphag_;
        }

        //- Return the phase-fraction of gaseous phase
        volScalarField& alphag()
        {
            return alphag_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
