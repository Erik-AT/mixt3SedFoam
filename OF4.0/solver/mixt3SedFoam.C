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

Application
    mixt3SedFoam

Description
	Three-phase model for sand transport in free surface flows
	Start:Wed Apr 25 17:10:09 2018

Author
    Mohamed Ouda

\*---------------------------------------------------------------------------*/
//General Finite Volume tools
#include "fvCFD.H"
#include "SortableList.H"

//For alphac equation
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"

//Transport, drag and turbulence models
#include "immiscibleIncompressibleThreePhaseMixture.H"
#include "mixture3DragModel.H"
#include "turbulenceModel.H"
#include "CompressibleTurbulenceModel.H"

//Numerical Schemes (check Y)
#include "gaussLaplacianScheme.H"
#include "uncorrectedSnGrad.H"
#include "localEulerDdtScheme.H"

//Solution algorithm 
#include "pimpleControl.H"

//Boundary conditions
#include "fixedValueFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"

#include "CorrectPhi.H"
#include "fvcSmooth.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "readGranularProperties.H"
    #include "correctPhi.H"

	// Validate the turbulence fields (mainly eddy viscosity) after construction of the turbulence 
	// model object from initial turbulence fields
    turbulence->validate();

    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "alphasCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
			//Read alpha & alphac solution controls            
			#include "alphasControls.H"

			//Solve the water and sediment volume conservation equations
            #include "alphasEqnSubCycle.H"

			//correct mixture viscosity and interface properties
            mixture.correct();

            //Calculate solid phase pressure and particle viscosity
            #include "granularRheologyCalc.H"		

			//Solve the momentum equation
            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
				//Solve the pressure equation
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
				//Solve turbulence closure equations and update the eddy viscosity
                turbulence->correct();
            }
        }

		//correct mixture viscosity
        mixture.correct();
		mum = mixture.mu();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
