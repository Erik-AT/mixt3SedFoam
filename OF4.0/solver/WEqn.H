{
	//Algebraic Slip Model (ASM) with turbulent dispersion and shear induced self-diffusion
	
	//Calculate the turbulent dispersion compoenent according to Burns et al.(2004)		
	if (turbulentDispersion)
	{

		Wtd = -turbulence->nut()/Sct	
									*(
										 fvc::grad(alphad)/max(alphad, scalar(1.0e-3) )
									   - fvc::grad(alphac)/max(alphac, scalar(1.0e-3) )
									 );

		//Keep only counter-gravity components
		Wtd*= 1.0 - pos(Wtd&g);
		Wtd.correctBoundaryConditions();

	}

	//Calculate the shear induced self-diffusion according to Zhang and Acrivos (1994)
	if (shearInduced)
	{
		#include "strainRates.H"
		Wsd = - 0.15*sqr(dd*alphad/max(gamma, scalar(1.0e-3)))*fvc::grad(magD);

		//Keep only counter-gravity components
		Wsd*= 1.0 - pos(Wsd&g);

		Wsd.correctBoundaryConditions();
	}

	//Scale the sediment pressure gradient to allow for larger time step
	//At the limit of alphad = alphaMax, the scaling factor = 1.0 
	volVectorField gradpsf(alphaMax*gradPsf/max(alphad, scalar(1.0e-3))); 

	// CD-W coupling loop
	for (int ite = 0; ite < nWcorr; ++ite)
	{
		Wi = 4.0/3.0*sqr(dd)
		*(rhod-rho)
		/(rhoc*mixture.nucModel().nu()*drag->CD()*drag->Re())
		*(g-((U)&fvc::grad(U))-fvc::ddt(U))();

		// if alphad < 1.0e-6 (one-way coupling), then set Wi = 0
		Wi*= pos(alphad-scalar(1.0e-6));

		Wi.correctBoundaryConditions();

		//update drag coefficient 
		drag->correct();
	}
	
	//Balance between particle pressure and particle weight 
	gradpNet = mag(gradpsf & gUnit) - alphad*(rhod-rhoc)*mag(g);
	Wi*= neg(gradpNet);	
	

	W = Wi;

    W.correctBoundaryConditions();

	W += (Wtd + Wsd);		
}


 




