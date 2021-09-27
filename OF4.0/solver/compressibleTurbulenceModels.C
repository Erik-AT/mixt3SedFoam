/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "CompressibleTurbulenceModel.H"
#include "immiscibleIncompressibleThreePhaseMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "makeTurbulenceModel.H"

#include "laminar.H"
#include "RASModel.H"
#include "LESModel.H"

makeBaseTurbulenceModel
(
    geometricOneField,
    volScalarField,
    compressibleTurbulenceModel,
    CompressibleTurbulenceModel,
    immiscibleIncompressibleThreePhaseMixture
);

#define makeRASModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (                                                                          \
        immiscibleIncompressibleThreePhaseMixtureCompressibleTurbulenceModel,   \
        RAS,                                                                   \
        Type                                                                   \
    )

#define makeLESModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (                                                                          \
        immiscibleIncompressibleThreePhaseMixtureCompressibleTurbulenceModel,   \
        LES,                                                                   \
        Type                                                                   \
    )


//Standard k-epsilon model
#include "kEpsilon.H"
makeRASModel(kEpsilon);

//Standard k-epsilon model + buoyancy term (OpenFOAM std library)
#include "buoyantKEpsilon.H"
makeRASModel(buoyantKEpsilon);

//Standard k-epsilon model + buoyancy term + drag term + (optional) free surface term
#include "mixtureKEpsilon.H"
makeRASModel(mixtureKEpsilon);

#include "twophaseMixingLength.H"
makeRASModel(twophaseMixingLength);


#include "Smagorinsky.H"
makeLESModel(Smagorinsky);

#include "kEqn.H"
makeLESModel(kEqn);

#include "dynamicKEqn.H"
makeLESModel(dynamicKEqn);


// ************************************************************************* //
