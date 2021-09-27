//Diffusion stresses
	word tauDmScheme("div(tauDm)"); 
	volSymmTensorField tauDm (
							  - alphad*alphac
							  * rhod*rhoc/max(gammaRhom, dimensionedScalar("small", dimDensity, 1.0e-3)) 
							  * sqr(W)       
							);

//Mixture effective viscosity
	volScalarField muEff (
							mixture.mu() //returns: gamma*mu_w + (1-gamma)mu_a
						  +	gamma * mup  //particulate viscosity
						  + mixture.rho() * turbulence->nut() //eddy viscosity
						 );

//Construct the Momentum equation 
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) 
	  + fvm::div(rhoPhi, U)
      - fvm::laplacian(muEff, U)
	  - fvc::div(muEff * dev2(T(fvc::grad(U))))
      - fvc::div(tauDm, tauDmScheme)
    );

    UEqn.relax();

	//Solve the momentum equation if required
    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
				  - fvc::snGrad(p_rgh)
                  - ghf*fvc::snGrad(rho)
                  + mixture.surfaceTensionForce()
                ) * mesh.magSf()
            )
            - gradPsf
            - gradPsk
        );

    }

