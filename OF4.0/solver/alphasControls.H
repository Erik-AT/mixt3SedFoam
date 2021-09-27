//Alphac Controls

	const dictionary& alphacControls = mesh.solverDict(alphac.name());

	label nAlphacCorr(readLabel(alphacControls.lookup("nAlphacCorr")));

	label nAlphaSubCycles(readLabel(alphacControls.lookup("nAlphaSubCycles")));

	bool MULESCorr(alphacControls.lookupOrDefault<Switch>("MULESCorr", false));

	// Apply the compression correction from the previous iteration
	// Improves efficiency for steady-simulations but can only be applied
	// once the alphac field is reasonably steady, i.e. fully developed
	bool alphacApplyPrevCorr
	(
		alphacControls.lookupOrDefault<Switch>("alphacApplyPrevCorr", false)
	);
	  
	// Isotropic compression coefficient
	scalar icAlpha
	(
		alphacControls.lookupOrDefault<scalar>("icAlpha", 0)
	);
	 
	// Shear compression coefficient
	scalar scAlpha
	(
		alphacControls.lookupOrDefault<scalar>("scAlpha", 0)
	);

//Alphad Controls

	const dictionary& alphadControls = mesh.solverDict(alphad.name());

	label nAlphadCorr(readLabel(alphadControls.lookup("nAlphadCorr")));
