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
    Foam::incompressibleThreePhaseMixture

Description
    A three-phase incompressible transportModel

SourceFiles
    incompressibleThreePhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleThreePhaseMixture_H
#define incompressibleThreePhaseMixture_H

#include "incompressible/transportModel/transportModel.H"
#include "incompressible/viscosityModels/viscosityModel/viscosityModel.H"
#include "mixture3ViscosityModel.H"
#include "threePhaseMixture.H"
#include "IOdictionary.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class incompressibleThreePhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class incompressibleThreePhaseMixture
:
    public IOdictionary,
    public transportModel,
    public threePhaseMixture
{
protected:

// Protected data


        //Pointer to suspension (dispersed + carrier) viscoity model	[new]
        autoPtr<mixture3ViscosityModel> muModel_;
       
        //Pointer to carrier phase viscoity model  						[nuModel1_]
        autoPtr<viscosityModel> nucModel_;
        
        //Pointer to gaseous phase viscosity model 						[nuModel2_]
        autoPtr<viscosityModel> nugModel_;
        

        //Density of dispersed phase									[new]	
        dimensionedScalar rhod_;
     
        //Density of carrier phase										[rho1_]
        dimensionedScalar rhoc_;
        
        //Density of gaseous phase										[rho2_]
        dimensionedScalar rhog_;
        

        //- Diameter of the dispersed phase particles					[new]
        dimensionedScalar dd_;

        //- Maximum dispersed phase-fraction (i.e. packing limit)		[new]
        scalar alphaMax_;
      
        //Turbulent Schmidt number										[new]
        dimensionedScalar Sct_;

		//Flag to apply packing limiter in aphad_Eqn					[new]
		Switch packingLimiter_; 

		//Maximum particle vicosity (for numerical stability)			[new]
		dimensionedScalar mupMax_;

		//Flag to apply limiter on the particle velocity				[new]
		Switch ULimiter_; 

        const volVectorField& U_;									//	[no change]
        const surfaceScalarField& phi_;								//	[no change]

        volScalarField nu_;											//	[no change]


    // Private Member Functions

        //- Calculate and return the laminar viscosity					[old]
        void calcNu();


public:

    TypeName("incompressibleThreePhaseMixture");


    // Constructors

        //- Construct from components
        incompressibleThreePhaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~incompressibleThreePhaseMixture()
    {}

	// Member Functions

        //- Return const-access to the suspension viscosityModel		[new]
        const mixture3ViscosityModel& muModel() const
        {
            return muModel_();
        }

        //- Return const-access to the carrier phase viscosityModel		[nuModel1()]
        const viscosityModel& nucModel() const
        {
            return nucModel_();
        }
        
        //- Return const-access to the gaseous phase viscosityModel		[nuModel2()]
        const viscosityModel& nugModel() const
        {
            return nugModel_();
        }

        //- Return const-access to the dispersed phase density			[new]
        const dimensionedScalar& rhod() const
        {
            return rhod_;
        }

        //- Return const-access to carrier phase density				[rho1()]				
        const dimensionedScalar& rhoc() const
        {
            return rhoc_;
        };
        
        //- Return const-access to gaseous phase density				[rho2()]
        const dimensionedScalar& rhog() const
        {
            return rhog_;
        };

        //- Return the diameter of the dispersed-phase particles		[new]
        const dimensionedScalar& dd() const
        {
            return dd_;
        }
        
        //- Return the turbulent Schmidt number							[new]
        const dimensionedScalar& Sct() const
        {
            return Sct_;
        }
        
		//- Return flag to apply packing limiter in aphad_Eqn			[new]
		const Switch& packingLimiter() const 
		{
			return packingLimiter_;
		}

		//- Return maximum particle vicosity (for numerical stability)	[new]
		const dimensionedScalar mupMax() const 
		{
			return mupMax_;
		}

		//- return flag to apply limiter on the particle velocity		[new]
		const Switch& ULimiter() const
		{
			return ULimiter_;
		} 

        //- Return the maximum phase fraction (i.e. packing limit)		[new]
        scalar alphaMax() const
        {
            return alphaMax_;
        }

        //- Return const-access to the mixture mass weighted velocity	[no change]
        const volVectorField& U() const
        {
            return U_;
        }



        //- Return the dynamic mixture laminar viscosity: suspension viscosity + gas viscosity.
        //  Dispersed phase viscosity is calculted separately using granular rheology
        tmp<volScalarField> mu() const;


		//- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> muf() const;

		//- Return the dynamic mixture viscosity for patch
        virtual tmp<scalarField> mu(const label patchi) const
        {
			volScalarField mu_(mu());
            return mu_.boundaryField()[patchi];
        }

        
        //- Return the kinematic mixture laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the kinematic mixture viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }   

		//- Return the face-interpolated kinematic laminar viscosity
        tmp<surfaceScalarField> nuf() const;     

        //- Return the mixture density based on limited alphac, to be used in mu & nu calculations
        virtual tmp<volScalarField> rhol() const
        {
			const volScalarField limitedAlphac
			(
				"limitedAlphac",
				min(max(alphac_, scalar(0)), scalar(1))
			);
            return alphad_*rhod_ + limitedAlphac*rhoc_ + (scalar(1) - alphad_ - limitedAlphac)*rhog_;
        }

        //- Return the mixture density
        virtual tmp<volScalarField> rho() const
        {
            return alphad_*rhod_ + alphac_*rhoc_ + (scalar(1) - alphad_ - alphac_)*rhog_;
        }


        //- Return the suspension density
 	    virtual tmp<volScalarField> rhom() const
		    {

            return (
						  (alphad_*rhod_ + alphac_*rhoc_)
						/ max(alphad_ + alphac_, scalar(1.0e-06))
						* neg(scalar(0) - alphad_ - alphac_)
					);
        }

       
        //- Return the suspension density * suspension volume fraction
        virtual tmp<volScalarField> gammaRhom() const
        {
            return alphad_*rhod_ + alphac_*rhoc_;
        }
        
	
        //- Correct the laminar viscosity: this function will corret both dynamic and kinematic 
		//- viscosity since calcNu() will call mu()
        virtual void correct()
        {

			calcNu();
        }

        //- Read base transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
