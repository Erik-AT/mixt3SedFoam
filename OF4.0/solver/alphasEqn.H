{
	//Numerical Schemes
	word alphadScheme("div(phi,alphad)");
    word alphadrScheme("div(phir,alphad)");
    word alphacScheme("div(phi,alphac)");
    word alphacrScheme("div(phirb,alphac)");

    // Standard face-flux compression coefficient
    surfaceScalarField phic(mixture.cAlpha()*mag(phi/mesh.magSf()));

    // Add the optional isotropic compression contribution
    if (icAlpha > 0)
    {
        phic *= (1.0 - icAlpha);
        phic += (mixture.cAlpha()*icAlpha)*fvc::interpolate(mag(U));
    }

    surfaceScalarField::Boundary& phicBf = phic.boundaryFieldRef();

    // Do not compress interface at non-coupled boundary faces
    // (inlets, outlets etc.)
    forAll(phic.boundaryField(), patchi)
    {
        fvsPatchScalarField& phicp = phicBf[patchi];

        if (!phicp.coupled())
        {
            phicp == 0;
        }
    }

	#include "ULimiter.H"  

	//Relative flux between dispersed and continuous phase(s)
	#include "WEqn.H"  

	surfaceScalarField phiW(fvc::flux(W));


	//### Alphac calaculations ####

	if(nAlphacCorr>0) //set nAlphacCorr= 0 for problems without free surface 
	{

		if (MULESCorr)
		{

			surfaceScalarField phi2
			(
				fvc::flux
				(
					phiW,
					- alphad / max(gamma, scalar(1e-03) )
					- alphad*alphac*(rhod-rhoc)/rho
					+ alphad*alphac*alphag*(1-S)/ max( gamma * (alphac + alphad*S),scalar(1e-03) )
				)
			);

			//Solve water mass conservation equation implicitly without the compression term
		    fvScalarMatrix alphacEqn
		    (
				 fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alphac)
		       + fv::gaussConvectionScheme<scalar>
		        (
		            mesh,
		            phi,
		            upwind<scalar>(mesh, phi)
		        ).fvmDiv(phi, alphac)
		       + fv::gaussConvectionScheme<scalar>
		        (
		            mesh,
		            phi2,
		            upwind<scalar>(mesh, phi2)
		        ).fvmDiv(phi2, alphac)
		    );

		    alphacEqn.solve();

		    Info<< "Carrier phase volume fraction = "
		        << alphac.weightedAverage(mesh.Vsc()).value()
		        << "  Min(" << alphac.name() << ") = " << min(alphac).value()
		        << "  Max(" << alphac.name() << ") = " << max(alphac).value()
		        << endl;

		    tmp<surfaceScalarField> talphacPhiUD(alphacEqn.flux());
		    alphacPhi = talphacPhiUD();

		    if (alphacApplyPrevCorr && talphacPhiCorr0.valid())
		    {
		        Info<< "Applying the previous iteration compression flux" << endl;
		        MULES::correct(alphac, alphacPhi, talphacPhiCorr0.ref(), 1, 0);

		        alphacPhi += talphacPhiCorr0();
		    }

		    // Cache the upwind-flux
		    talphacPhiCorr0 = talphacPhiUD;

			gamma = min(max(alphad + alphac, scalar(0)), scalar(1.0));
		
			alphac = gamma - alphad;  

		    alphag = 1 - gamma; 

			//Replace air in water-sand interface by water
			if(addWater)
			{
				alphac += alphag*(pos(alphad-alphag)*pos(alphac-0.1)*pos(0.99-alphac));
			}

			//Remove air from water-sand interface
			if(removeAir)
			{
				alphag *= (1 - pos(alphad-alphag)*pos(alphac-0.1)*pos(0.99-alphac) );
			}

			mixture.correct();
		}


		for (int aCorr=0; aCorr<nAlphacCorr; aCorr++)
		{
		    surfaceScalarField phir(phic*mixture.nHatf());

			tmp<surfaceScalarField> talphacPhiUn
		    (
		        fvc::flux
		        (
		            phi,
		            alphac,
		            alphacScheme
		        )
		      + fvc::flux
		        (
		           -fvc::flux( -phir, alphag, alphacrScheme),
		            alphac,
		            alphacrScheme
		        )
		      + fvc::flux
		        (
					   -fvc::flux( -phiW,
									- alphad / max(gamma, scalar(1e-03))
									- alphad*alphac*(rhod-rhoc)/rho
									+ alphad*alphac*alphag*(1-S)/ max( gamma * (alphac + alphad*S),scalar(1e-03) ),
									alphacrScheme
								 ),
					alphac,
					alphacrScheme
				)   
		    );


		    if (MULESCorr)
		    {
		        tmp<surfaceScalarField> talphacPhiCorr(talphacPhiUn() - alphacPhi);
		        volScalarField alphac0("alphac0", alphac);

		        MULES::correct(alphac, talphacPhiUn(), talphacPhiCorr.ref(), 1, 0);

		        // Under-relax the correction for all but the 1st corrector
		        if (aCorr == 0)
		        {
		            alphacPhi += talphacPhiCorr();
		        }
		        else
		        {
		            alphac = 0.5*alphac + 0.5*alphac0;
		            alphacPhi += 0.5*talphacPhiCorr();
		        }
		    }
		    else
		    {
		        alphacPhi = talphacPhiUn;

		        MULES::explicitSolve(alphac, phi, alphacPhi, 1, 0);
		    }

			gamma = min(max(alphad + alphac, scalar(0)), scalar(1.0));

			alphac = gamma - alphad; 

		    alphag = 1 - gamma; 

			//Replace air in water-sand interface by water
			if(addWater)
			{
				alphac += alphag*(pos(alphad-alphag)*pos(alphac-0.1)*pos(0.99-alphac));
			}

			//Remove air from water-sand interface
			if(removeAir)
			{
				alphag *= (1 - pos(alphad-alphag)*pos(alphac-0.1)*pos(0.99-alphac) );
			}

			mixture.correct();
		}

		if (alphacApplyPrevCorr && MULESCorr)
		{
		    talphacPhiCorr0 = alphacPhi - talphacPhiCorr0;
		}

		Info<< "Carrier phase volume fraction = "
		    << alphac.weightedAverage(mesh.Vsc()).value()
		    << "  Min(" << alphac.name() << ") = " << min(alphac).value()
		    << "  Max(" << alphac.name() << ") = " << max(alphac).value()
		    << endl;
	}
	
	//For problems without free surface: 
	else
	{
		Info<< "Free surface calculation is turned OFF" << endl << endl; 
		gamma  = 1.0;
		alphac = 1.0 - alphad;
		alphag = 0;
	}


	//### Alphad calculations ###
    for (int acorr=0; acorr<nAlphadCorr; acorr++)
    {

		surfaceScalarField phi3 
		(
			 -fvc::flux( -phiW,
							   alphac / max(gamma, scalar(1e-03) )
						 	 - alphad*alphac*(rhod-rhoc)/rho
							 + alphad*alphac*alphag*(1-S)/ max( gamma * (alphac + alphad*S),scalar(1e-03) ), 
							 alphadrScheme
					   )
		);


		//Fully implicit 
        fvScalarMatrix alphadEqn
        (
             fvm::ddt(alphad)
           + fvm::div(phi , alphad, alphadScheme)
		   + fvm::div(phi3, alphad, alphadrScheme)
        );

	    alphadEqn.relax();
	    alphadEqn.solve();
	

		//Limit alphad to its maximum value without violating the sediment mass conservation
	    if (mixture.packingLimiter())
	    {
		int i=0;
	        while (i<10 && max(alphad).value()>=alphaMax-0.02)
	        {
	            #include "phaseFractionLimiter.H"
	            i++;
	        }
	    }
		

		gamma = min(max(alphad + alphac, scalar(0)), scalar(1.0));

		alphac = gamma - alphad;  

        alphag = 1 - gamma; 

		//Replace air in water-sand interface by water
		if(addWater)
		{
			alphac += alphag*(pos(alphad-alphag)*pos(alphac-0.1)*pos(0.99-alphac));
		}

		//Remove air from water-sand interface
		if(removeAir)
		{
			alphag *= (1 - pos(alphad-alphag)*pos(alphac-0.1)*pos(0.99-alphac) );
		}

		alphadPhi = alphadEqn.flux();

		//Update mixture density and viscosity		
		#include "updateRho.H"
		mixture.correct();

        Info<< "Dispersed Phase volume fraction = "
            << alphad.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alphad.name() << ") = " << min(alphad).value()
            << "  Max(" << alphad.name() << ") = " << max(alphad).value()
            << endl;
	}

    // Calculate the end-of-time-step mass flux
    rhoPhi = alphadPhi*(rhod - rhog) + alphacPhi*(rhoc - rhog) + phi*rhog;

	//Update alphas field
	alphas = 0.0*alphac + 1.0*alphag + 2.0*alphad;

	//Check for sediment and water mass conservation
	sedVol = gSum((mesh.V()*alphad.internalField())());
	waterVol = gSum((mesh.V()*alphac.internalField())());

	alphaMaxActual = (max(alphad.internalField())).value();
	alphaMinActual = (min(alphad.internalField())).value();

	std::ofstream file;
	file.open ("MassConservation.txt", std::ofstream::out | std::ofstream::app);
	file << runTime.timeName() << "	" << sedVol << "	" << alphaMaxActual << "	" 
		 << alphaMinActual << "	" << waterVol << std::endl;
	file.close();


}
