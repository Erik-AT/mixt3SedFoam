{
	//****************************** Particle pressure calculations *******************************           

	// A- Calculate the enduring contact particle component of the pressure (psf)
	// according to Johnson and Jackson (1987)
	psf =
	(
		Cp*pow(max(alphad - alphaMinFriction, scalar(0)), C1)
		/pow(max(alphaMax - alphad, scalar(1.0e-20)), C2)
	);
	
	//Update total particle pressure
	ps = psf;


	// Gradient of psf: Dealing psf as body force or source term and calculate 
	// its gradient in two steps in line with OpenFOAM standards:

		//1- Calcualte the the face flux of gradpsf
		gradPsff = fvc::snGrad(psf)*mesh.magSf();

		//If pressure B.C. is set to zero gradient, then set gradpsff to zero	
		forAll(p.boundaryField(), patchi)
		{
			if (isA<zeroGradientFvPatchScalarField>(p.boundaryField()[patchi]))
			{
				gradPsff.boundaryFieldRef()[patchi] = 0.0;
			}
		}

		//2- Reconstruct gradpsf at cell center from its face flux (gradPsff)
		gradPsf = (fvc::reconstruct(gradPsff));


	
	// B- Calculate the collisional/kinetic component of the particle pressure (psk) 
	

	#include "strainRates.H"

	//Calculate Stokes number (Cassar et al., 2005)
	volScalarField St(dd*sqrt(rhod*psk)/(rhoc*mixture.nucModel().nu()));

/*
	//Calculate Stokes number (unified approach)
	volScalarField St(sqr(dd)*rhod*magD/(rhoc*mixture.nucModel().nu())); 
*/

	if (kineticPressure)
	{
		//Store old value of psk
		const volScalarField pskOld(psk);
			
		//Calculate new value for psk (Chauchat, 2017)
		psk = ( pos(St-1)* rhod*sqr(dd)* magsqrD + neg(St-1)*rhoc*mixture.nucModel().nu()* magD )
			  * sqr(Bphi*alphad/max(alphaMax - alphad, scalar(1.0e-3)));

/*
		//Calculate new value for psk (unified approach)
		psk = rhoc*mixture.nucModel().nu()*(1+St)*sqr(Bphi*alphad/max(alphaMax - alphad, scalar(1.0e-3)))*magD;
*/

		//Relax psk 
		psk = pskOld + relaxPsk*(psk-pskOld);
	
		//Face flux of gradpsk 		
		 gradPskf = fvc::snGrad(psk)*mesh.magSf();
	   
		forAll(p.boundaryField(), patchi)
		{
			if (isA<zeroGradientFvPatchScalarField>(p.boundaryField()[patchi]))
			{
				gradPskf.boundaryFieldRef()[patchi] = 0.0;
			}
		}

		//Gradient of psk
		gradPsk = (fvc::reconstruct(gradPskf));
	
		//Add psk to the total particle pressure
		ps+= psk;		
	
	}
	
	//*************************** Particle shear viscosity calculations ***************************     

	//Calculate Mu(I) based on the flow regime (Chauchat, 2017)
	volScalarField muI 
					(	
						mus + (mu2 - mus)*magD/
							(
								pos(St-1)*I0*sqrt(ps/rhod)/dd 	
							  +	neg(St-1)*I0*ps/(rhoc*mixture.nucModel().nu())
							  + magD + lambdaR
							)
					);

/*
	//Calculate Mu(I) (unified approach)
	volScalarField muI 
					(	
						mus + (mu2 - mus)*magD/
							( 	
							  I0*ps/( (1+St)*rhoc*mixture.nucModel().nu() ) + magD + lambdaR
							)
					);
*/
		  
	//Calculate the particle shear viscosity
	mup = muI*ps/sqrt(magsqrD + sqr(lambdaR));

	//Calculate the particle viscosity on boundaries
    forAll(alphad.boundaryField(), patchi)
    {
        if (isA<zeroGradientFvPatchScalarField>(alphad.boundaryField()[patchi]))
        {
            mup.boundaryFieldRef()[patchi] =
            (
                muI.boundaryField()[patchi]*ps.boundaryField()[patchi]
               /sqrt(magsqrD.boundaryField()[patchi] + sqr(lambdaR).value())
            );
        }
    }

	//Limit the particle viscosity for numerical stability
	mup.min(mupMax);
}
