
const dictionary& granularDict(transProperties.subDict("granularRheology"));

//Parameters for frictional pressure of Johnson and Jackson (1987)
 
	scalar alphaMax(mixture.alphaMax());

	dimensionedScalar alphaMinFriction
	(
		granularDict.lookupOrDefault
		(
			"alphaMinFriction",
			dimensionedScalar("alphaMinFriction",dimensionSet(0,0,0,0,0,0,0), 0.523598 )
		)
	);

	dimensionedScalar Cp
	(
		granularDict.lookupOrDefault
		(
			"Cp",
			dimensionedScalar("Cp",dimensionSet(1,-1,-2,0,0,0,0), 0.05 )
		)
	);


	dimensionedScalar C1
	(
		granularDict.lookupOrDefault
		(
			"C1",
			dimensionedScalar("C1",dimensionSet(0,0,0,0,0,0,0), 3 )
		)
	);

	dimensionedScalar C2
	(
		granularDict.lookupOrDefault
		(
			"C2",
			dimensionedScalar("C2",dimensionSet(0,0,0,0,0,0,0), 5 )
		)
	);  

//Parameters for kinetic pressure
	bool  kineticPressure (granularDict.lookupOrDefault("kineticPressure", true));
	dimensionedScalar Bphi(granularDict.lookup("Bphi"));
	dimensionedScalar relaxPsk(granularDict.lookup("relaxPsk"));
	
//Parameters for particle shear viscosity
	dimensionedScalar lambdaR
	(
		granularDict.lookupOrDefault
		(
			"lambdaR",
			dimensionedScalar("lambdaR",dimensionSet(0,0,-1,0,0,0,0), 1.0e-6 )
		)
	);
	
	dimensionedScalar mus(granularDict.lookup("mus"));
	
	dimensionedScalar mu2(granularDict.lookup("mu2"));
	
	dimensionedScalar I0(granularDict.lookup("I0"));


