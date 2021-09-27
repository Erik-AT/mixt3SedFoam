//Calculate shear rate using mixture velocity
	volSymmTensorField D = dev(symm(fvc::grad(U)));
	volScalarField magD = ::sqrt(2.0)*mag(D);
	volScalarField magsqrD = sqr(magD);

