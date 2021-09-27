//Read transport properties dictionary
IOdictionary transProperties
(
	IOobject
	(
		"transportProperties",
		runTime.constant(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	)
);


//Creating driving pressure gradient
	Info<< "\nReading transportProperties for gradP\n" <<endl;

	//Oscillating pressure gradient
	dimensionedVector gradPOSC("gradPOSC", dimensionSet(1, -2, -2, 0, 0, 0, 0), vector(0, 0, 0));

	//Mean pressure gradient
	dimensionedVector gradPMean
	(
		transProperties.lookupOrDefault
		(
			"gradPMean",
			dimensionedVector("gradPMean", dimensionSet(1, -2, -2, 0, 0, 0, 0), vector(0, 0, 0)) 
		)
	);

	// Pressure oscillation amplitude 1
	dimensionedVector gradPAMP1
	(
		transProperties.lookupOrDefault
		(
			"gradPAMP1",
			dimensionedVector("gradPAMP1", dimensionSet(1, -2, -2, 0, 0, 0, 0), vector(0, 0, 0)) 
		)
	);

	// Pressure oscillation amplitude 2 (in case of uisng Stokes 2nd order)
	dimensionedVector gradPAMP2
	(
		transProperties.lookupOrDefault
		(
			"gradPAMP2",
			dimensionedVector("gradPAMP2", dimensionSet(1, -2, -2, 0, 0, 0, 0), vector(0, 0, 0)) 
		)
	);

	// Pressure oscillation period
	dimensionedScalar oscpT
	(
		transProperties.lookupOrDefault
		(
			"oscpT",
			dimensionedScalar("oscpT", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1.0) 
		)
	);

	//Initial deg (phase shift)
	dimensionedScalar initTheta
	(
		transProperties.lookupOrDefault
		(
			"initTheta",
			dimensionedScalar("initTheta", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0) 
		)
	);



//Names of different phases
word dispersedPhaseName(wordList(transProperties.lookup("phases"))[0]);
word carrierPhaseName(wordList(transProperties.lookup("phases"))[1]);
word gasPhaseName(wordList(transProperties.lookup("phases"))[2]);


//Modified pressure = static pressure + hydrostatic pressure (body force)
Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//Frictional or enduring contact component of the particle pressure
volScalarField psf
(
    IOobject
    (
        "psf",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT, //Allow for initialization of Psf
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(1 , -1 ,-2 ,0 , 0), 0.0)
);



// Gradient of psf (will be used in the momentum equation)
volVectorField gradPsf
(
    IOobject
    (
        "gradPsf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(1, -2, -2, 0, 0), vector::zero)
);

//Face flux of gradpsf (will be used in the pressure equation)
surfaceScalarField gradPsff(fvc::snGrad(psf)*mesh.magSf());

//Kinetic or collisional component of the particle pressure
Info<< "Reading field psk\n" << endl;
volScalarField psk
(
    IOobject
    (
        "psk",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//Gradient of psk (will be used in the momentum equation)
volVectorField gradPsk
(
    IOobject
    (
        "gradPsk",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(1, -2, -2, 0, 0), vector::zero)
);

//Face flux of gradpsk (will be used in the pressure equation)
surfaceScalarField gradPskf(fvc::snGrad(psk)*mesh.magSf()); 

//Total particle pressure
volScalarField ps(psf + psk);

//The particle viscosity
volScalarField mup
(
    IOobject
    (
        "mup",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(1 , -1 ,-1 ,0 , 0), 0.0)
);


//Mass averaged mixture velocity 
Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//phi = U_f.S_f; the face volumetric flux based on U
#include "createPhi.H"

//Mixture transport model
Info<< "Reading immiscibleIncompressibleThreePhaseMixture\n" << endl;
immiscibleIncompressibleThreePhaseMixture mixture(U, phi);

//Volume fraction of:
	//- Dispersed phase
	volScalarField& alphad(mixture.alphad());
	//- Suspension
	volScalarField& gamma(mixture.gamma());
	gamma.write();
	//- Carrier phase
	volScalarField& alphac(mixture.alphac());
	//- Gaseous phase
	volScalarField& alphag(mixture.alphag());
	alphag.write();

//Alphas field for visualization purpose
volScalarField alphas
(
    IOobject
    (
        "alphas",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	  0.0*alphac + 1.0*alphag + 2.0*alphad 
);
alphas.write();

//Density of:
	//- Dispersed phase	
	const dimensionedScalar& rhod = mixture.rhod();
	//- Carrier phase	
	const dimensionedScalar& rhoc = mixture.rhoc();
	//- Gaseous phase	
	const dimensionedScalar& rhog = mixture.rhog();
	
	const dimensionedScalar S = rhod/rhoc;


//Dispersed phase diameter
const dimensionedScalar& dd(mixture.dd());

//Turbulent Schmidt number
const dimensionedScalar& Sct(mixture.Sct());

//Maximum particle viscosity
dimensionedScalar mupMax( mixture.mupMax() );

// Mixture density
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mixture.rho()
);

// Need to store rho for ddt(rho, U)
rho.oldTime();

//Suspension density
volScalarField rhom
(
    IOobject
    (
        "rhom",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mixture.rhom()
);

//gammaRhom
volScalarField gammaRhom
(
    IOobject
    (
        "gammaRhom",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mixture.gammaRhom()
);

//Suspension viscosity
volScalarField mum
(
    IOobject
    (
        "mum",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mixture.mu()
);

// Mass flux = rho*phi
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);


//Relative velocity between water and sediment: W = v_s - v_w
#include "createW.H"

//Components of W
volVectorField Wi
(
    IOobject
    (
        "Wi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    W*0
);

volVectorField Wtd
(
    IOobject
    (
        "Wtd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    W*0
);

volVectorField Wsd
(
    IOobject
    (
        "Wsd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    W*0
);

// Construct compressible turbulence model
autoPtr <CompressibleTurbulenceModel <immiscibleIncompressibleThreePhaseMixture> > turbulence
(
    CompressibleTurbulenceModel<immiscibleIncompressibleThreePhaseMixture>
    ::New(rho, U, rhoPhi, mixture)
);

//Construct Drag model
autoPtr<mixture3DragModel> drag(mixture3DragModel::New(mixture, Wi));

//mixture3DragModel& dragModel = dragModelPtr();


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

//unit vector in the gravity direction		
vector gUnit = (g / mag(g)).value();

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

//The gradpsf - alphad*(rhod-rhoc)*g
volScalarField gradpNet
(
    IOobject
    (
        "gradpNet",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(1 , -2 ,-2 ,0 , 0), 0.0)
);


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alphad.name());
mesh.setFluxRequired(alphac.name());

// MULES flux from previous time-step
surfaceScalarField alphacPhi
(
    IOobject
    (
        "alphacPhi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    phi*fvc::interpolate(alphac)
);

surfaceScalarField alphadPhi
(
    IOobject
    (
        "alphadPhi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    phi*fvc::interpolate(alphad)
);

// MULES Correction
tmp<surfaceScalarField> talphacPhiCorr0;

//The following two flags are used to solve the stability problem which occurs when air is entrained
//on water-sediment interface. It may cause mass conservation error which can be large in some cases 

	//Flag to remove air from water-sand interface
	Switch removeAir(transProperties.lookupOrDefault("removeAir", false));

	//Flag to add water instead of the removed air
	Switch addWater(transProperties.lookupOrDefault("addWater", false)); 

//Flag for stationary bed: 0 = stationary region
volScalarField packedInd
(
	IOobject
	(
	    "packedInd",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("One", dimless, 1.0)
);

//scalar alphaLimit(readScalar(transProperties.lookup("alphaLimit")));

//Check for sediment and water mass conservation
scalar sedVol = 0;
scalar waterVol = 0;
scalar alphaMaxActual = 0;
scalar alphaMinActual = 0;

std::ofstream file;
file.open ("MassConservation.txt", std::ofstream::out | std::ofstream::app);
file << "Time" << "	" << "sedVol" << "	" << "alphaMax" << "	" << "alphaMin" << "	" << "waterVol"<< std::endl;


const labelUList & Owner(mesh.owner());
const labelUList & Neigh(mesh.neighbour());
