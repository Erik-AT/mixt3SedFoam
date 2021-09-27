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

Global
    alphasCourantNo

Description
    Description
    Calculates and outputs the mean and maximum interface Courant Numbers for:
	1- disperesed-carrier interface
	2- carrier-gaseous interface


\*---------------------------------------------------------------------------*/

scalar maxAlphaCo
(
    readScalar(runTime.controlDict().lookup("maxAlphaCo"))
);

scalar alphadCoNum = 0.0;
scalar alphacCoNum = 0.0;

scalar meanAlphadCoNum = 0.0;
scalar meanAlphacCoNum = 0.0;

if (mesh.nInternalFaces())
{
     scalarField sumPhi1
     (
         pos(alphad - 0.01)*pos(0.99 - alphad)			//disperesed-carrier interface
        *fvc::surfaceSum(mag(phi))().primitiveField()
     );

    scalarField sumPhi2
    (
        mixture.nearInterface()().primitiveField()		//carrier-gaseous interface
       *fvc::surfaceSum(mag(phi))().primitiveField()
    );

	alphadCoNum = 0.5*gMax(sumPhi1/mesh.V().field())*runTime.deltaTValue();
	alphacCoNum = 0.5*gMax(sumPhi2/mesh.V().field())*runTime.deltaTValue();

     meanAlphadCoNum =
         0.5*(gSum(sumPhi1)/gSum(mesh.V().field()))*runTime.deltaTValue();
     meanAlphacCoNum =
        0.5*(gSum(sumPhi2)/gSum(mesh.V().field()))*runTime.deltaTValue();
}

Info<< dispersedPhaseName << "-" << carrierPhaseName << " Interface Courant Number mean: " << meanAlphadCoNum
    << " max: " << alphadCoNum << endl;

Info<< carrierPhaseName << "-" << gasPhaseName << " Interface Courant Number mean: " << meanAlphacCoNum
    << " max: " << alphacCoNum << endl << endl;

// ************************************************************************* //
