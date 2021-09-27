/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "incompressibleThreePhaseMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleThreePhaseMixture, 0);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::incompressibleThreePhaseMixture::calcNu()
{
    nucModel_->correct();
    nugModel_->correct();

   /* const volScalarField limitedGamma
    (
        "limitedGamma",
        min(max(gamma_, scalar(0)), scalar(1))
    );*/

    // Average kinematic viscosity calculated from dynamic viscosity
    nu_ = mu()/rhol();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleThreePhaseMixture::incompressibleThreePhaseMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    threePhaseMixture(U.mesh(), *this),

    muModel_
    (
        mixture3ViscosityModel::New
        (
            "mu",
            subDict(phase1Name_),
            U
        )
    ),

    nucModel_
    (
        viscosityModel::New
        (
            "nuc",
            subDict(phase2Name_),
            U,
            phi
        )
    ),
    
    nugModel_
    (
        viscosityModel::New
        (
            "nug",
            subDict(phase3Name_),
            U,
            phi
        )
    ),    


    rhod_("rho", dimDensity, muModel_->viscosityProperties()),
    rhoc_("rho", dimDensity, nucModel_->viscosityProperties()),
    rhog_("rho", dimDensity, nugModel_->viscosityProperties()),


    dd_
    (
        "d",
        dimLength,
        muModel_->viscosityProperties().lookup("d")
    ),
    
    alphaMax_(readScalar(muModel_->viscosityProperties().lookup("alphaMax"))),
    
    Sct_
    (
		"Sct", 
		dimless, 
		muModel_->viscosityProperties().lookupOrDefault("Sct", 0.7)
	),

	packingLimiter_
	(
		muModel_->viscosityProperties().lookup("packingLimiter")
	),

	mupMax_
	(
		muModel_->viscosityProperties().lookup("mupMax")
	),

	ULimiter_
	(
		muModel_->viscosityProperties().lookup("ULimiter")
	),

    U_(U),

    phi_(phi),

    nu_
    (
        IOobject
        (
            "nu",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedScalar("nu", dimViscosity, 0),
        calculatedFvPatchScalarField::typeName
    )

{
    calcNu();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::incompressibleThreePhaseMixture::mu() const
{
    const volScalarField limitedAlphac
    (
        min(max(alphac_, scalar(0)), scalar(1))
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "mu",
            (alphad_ + limitedAlphac)*rhoc_*muModel_->nu( nucModel_->nu(), alphaMax_ ) 
          + (scalar(1) - alphad_ - limitedAlphac)*rhog_*nugModel_->nu()
        )
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleThreePhaseMixture::muf() const
{
    const surfaceScalarField alphacf
    (
        min(max(fvc::interpolate(alphac_), scalar(0)), scalar(1))
    );

	const surfaceScalarField alphadf(fvc::interpolate(alphad_));

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "muf",
            (alphadf + alphacf)*rhoc_*fvc::interpolate( muModel_->nu(nucModel_->nu(), alphaMax_) )
          + (scalar(1) - alphadf - alphacf)*rhog_*fvc::interpolate(nugModel_->nu())
        )
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleThreePhaseMixture::nuf() const
{
    const surfaceScalarField alphacf
    (
        min(max(fvc::interpolate(alphac_), scalar(0)), scalar(1))
    );

	const surfaceScalarField alphadf(fvc::interpolate(alphad_));

//	const surfaceScalarField rhof (fvc::interpolate(rho()));

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "nuf",
            (
            (alphadf + alphacf)*rhoc_*fvc::interpolate( muModel_->nu(nucModel_->nu(), alphaMax_) )
          + (scalar(1) - alphadf - alphacf)*rhog_*fvc::interpolate(nugModel_->nu())
            )/(alphadf*rhod_ + alphacf*rhoc_ + (scalar(1) - alphadf - alphacf)*rhog_) 
        )
    );
}

bool Foam::incompressibleThreePhaseMixture::read()
{
    if (regIOobject::read())
    {
        if
        (
            muModel_().read(subDict(phase1Name_))
         && nucModel_().read(subDict(phase2Name_))
         && nugModel_().read(subDict(phase3Name_))
        )
        {
            muModel_->viscosityProperties().lookup("rho") >> rhod_;
            nucModel_->viscosityProperties().lookup("rho") >> rhoc_;
            nugModel_->viscosityProperties().lookup("rho") >> rhog_;
            

            dd_ = dimensionedScalar
            (
                "d",
                dimLength,
                muModel_->viscosityProperties().lookup("d")
            );

            alphaMax_ = readScalar(muModel_->viscosityProperties().lookup("alphaMax"));

			packingLimiter_ = Switch( muModel_->viscosityProperties().lookup("packingLimiter"));

			mupMax_ = muModel_->viscosityProperties().lookup("mupMax");

			ULimiter_ = Switch(muModel_->viscosityProperties().lookup("ULimiter"));

            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
