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

#include "mixtureKEpsilon.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void mixtureKEpsilon<BasicTurbulenceModel>::correctNut()
{
/*	//Limit nut at high strain zones according to P. A. Durbin (2009)
	
	dimensionedScalar Dsmall("lambdaR",dimensionSet(0,0,-1,0,0,0,0), 1.0e-6 );
    const volVectorField& U = this->U_;
	tmp<volScalarField> magD(::sqrt(2.0)*mag(dev(symm(fvc::grad(U)))));
    this->nut_ = min(Cmu_*sqr(k_)/epsilon_, ::sqrt(2.0/3.0)*k_/(magD() + Dsmall));
*/
	this->nut_ = Cmu_*sqr(k_)/epsilon_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}



template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> mixtureKEpsilon<BasicTurbulenceModel>::kSource() const
{
    const uniformDimensionedVectorField& g =
        this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");

    if (mag(g.value()) > SMALL)
    {
        return -fvm::SuSp(Gcoef() + FScoef() + Dcoef(), k_);
    }

    else
    {
		return tmp<fvScalarMatrix>
		(
		    new fvScalarMatrix
		    (
		        k_,
		        dimVolume*this->rho_.dimensions()*k_.dimensions()
		        /dimTime
		    )
    	);
    }

}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> mixtureKEpsilon<BasicTurbulenceModel>::epsilonSource() const
{
    const uniformDimensionedVectorField& g =
        this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");

    if (mag(g.value()) > SMALL)
    {
        vector gHat(g.value()/mag(g.value()));

        volScalarField v(gHat & this->U_);
        volScalarField u
        (
            mag(this->U_ - gHat*v)
          + dimensionedScalar("SMALL", dimVelocity, SMALL)
        );

        return -fvm::SuSp(C1_*tanh(mag(v)/u)*Gcoef() + FScoef() + C3eps_*Dcoef(), epsilon_);
    }

    else
    {
		return tmp<fvScalarMatrix>
		(
		    new fvScalarMatrix
		    (
		        epsilon_,
		        dimVolume*this->rho_.dimensions()*epsilon_.dimensions()
		        /dimTime
		    )
		);
    }

}

//Buoyancy damping term according to Henkes et al. (1991)
template<class BasicTurbulenceModel>
tmp<volScalarField>
mixtureKEpsilon<BasicTurbulenceModel>::Gcoef() const
{
    const uniformDimensionedVectorField& g =
        this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");

	const volScalarField& gamma(this->mesh_.objectRegistry::template
								 lookupObject<volScalarField>("gamma"));


    return
        (Cg_*Cmu_)*gamma*k_*(g & fvc::grad(this->rho_))*pos((g & fvc::grad(this->rho_)))
       /(epsilon_ + this->epsilonMin_);

}


//Drag damping term according to Hsu et al. (2004)
template<class BasicTurbulenceModel>
tmp<volScalarField>
mixtureKEpsilon<BasicTurbulenceModel>::Dcoef() const
{

	const volScalarField& alphad(this->mesh_.objectRegistry::template
								 lookupObject<volScalarField>(alphaName_));

	const volScalarField& Cd (this->mesh_.objectRegistry::template
							  lookupObject<volScalarField>("CD"));

	const volVectorField& Wi (this->mesh_.objectRegistry::template
							  lookupObject<volVectorField>("Wi"));

	//drag parameter
	volScalarField KK(0.75*Cd*rhoc_*mag(Wi)/dd_);
	
	//particle time scale
	volScalarField tp( rhod_/max(KK*(1-alphad), dimensionedScalar("small", dimensionSet(1,-3,-1,0,0,0,0), 1.0e-3)));

	//fluid time scale 
	volScalarField tf(k_/(6*(epsilon_ + this->epsilonMin_)));

	//Stokes number
	volScalarField St(tp/max(tf, dimensionedScalar("small", dimensionSet(0,0,1,0,0,0,0), 1.0e-6)));

	//correlation parameter
	volScalarField tmf(Foam::exp(-B_*St));

	return (2.0*KK*(1-tmf)*alphad); 
}


//Free surface damping term according to Naot and Rodi (1982)
template<class BasicTurbulenceModel>
tmp<volScalarField>
mixtureKEpsilon<BasicTurbulenceModel>::FScoef() const
{

	if (dampFS_)
	{
		const volScalarField& gamma(this->mesh_.objectRegistry::template
									 lookupObject<volScalarField>("gamma"));
		volScalarField limitedGamma = min(max(gamma, scalar(0)), scalar(1));
		const volVectorField gradGamma(fvc::grad(limitedGamma));

		return
		    2.5*pow(Cmu_,-0.25)*pow(k_,0.5)*this->rho_/(0.07*h_)
			*pos(mag(gradGamma));
	}

    else
    {
		return tmp<volScalarField>
		(
		    new volScalarField
		    (
		    IOobject
			(
			    "Zero",
			    this->runTime_.timeName(),
			    this->mesh_,
			    IOobject::NO_READ,
			    IOobject::NO_WRITE
			),
		    this->mesh_,
			dimensionedScalar("Zero", dimDensity/dimTime, 0)
		    )
		);
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
mixtureKEpsilon<BasicTurbulenceModel>::mixtureKEpsilon
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            1.92
        )
    ),
    C3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3",
            this->coeffDict_,
            -0.33
        )
    ),
    sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.3
        )
    ),

    Cg_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cg",
            this->coeffDict_,
            1.0
        )
    ),

    B_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "B",
            this->coeffDict_,
            0.25
        )
    ),

    C3eps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3eps",
            this->coeffDict_,
            1.2
        )
    ),

    C3k_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3k",
            this->coeffDict_,
            1.0
        )
    ),

	transProperties_
	(
		IOobject
		(
			"transportProperties",
			U.time().constant(),
			U.db(),
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	),

	dispersedName_(wordList(transProperties_.lookup("phases"))[0]),

	carrierPhaseName_(wordList(transProperties_.lookup("phases"))[1]),
	
	alphaName_(IOobject::groupName("alpha", dispersedName_)),

	rhoc_("rho", dimDensity, transProperties_.subDict(carrierPhaseName_)),
	rhod_("rho", dimDensity, transProperties_.subDict(dispersedName_)),
	dd_("d", dimLength, transProperties_.subDict(dispersedName_)),

    dampFS_
    (
        Switch::lookupOrAddToDict
        (
            "dampFS",
            this->coeffDict_,
            false
        )
    ),

	h_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "h",
            this->coeffDict_,
            dimLength,
			1.0
        )
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool mixtureKEpsilon<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());
		Cg_.readIfPresent(this->coeffDict());
		B_.readIfPresent(this->coeffDict());
		C3eps_.readIfPresent(this->coeffDict());
		C3k_.readIfPresent(this->coeffDict());


        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void mixtureKEpsilon<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;


    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))().v()
    );

    tmp<volTensorField> tgradU = fvc::grad(U);

	//turbulence shear production term
    volScalarField::Internal G
    (
        this->GName(),
        nut.v()*(dev(twoSymm(tgradU().v())) && tgradU().v())
    );
    tgradU.clear();

    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        C1_*alpha()*rho()*G*epsilon_()/k_()
      - fvm::SuSp(((2.0/3.0)*C1_ + C3_)*alpha()*rho()*divU, epsilon_)
      - fvm::Sp(C2_*alpha()*rho()*epsilon_()/k_(), epsilon_)
      + epsilonSource()
    );

    epsEqn.ref().relax();
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    bound(epsilon_, this->epsilonMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        C3k_*alpha()*rho()*G  //this term was modified
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*epsilon_()/k_(), k_)
      + kSource()
    );

    kEqn.ref().relax();

    solve(kEqn);

    bound(k_, this->kMin_);

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
