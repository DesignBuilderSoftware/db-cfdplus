/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Derived from OpenFOAM's MachNo function object.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2013-2018 FSD blueCAPE Lda  http://www.bluecape.com.pt
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "temperatureCelcius.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(temperatureCelcius, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        temperatureCelcius,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::temperatureCelcius::calc()
{
    if
    (
        foundObject<volScalarField>(fieldName_)
    )
    {
        const volScalarField& T = lookupObject<volScalarField>(fieldName_);

        Info<< "Minimum temperature: " << gMin(T) << endl;
        Info<< "Maximum temperature: " << gMax(T) << endl;

        return store
        (
            resultName_,
            T-zero_celsius_in_kelvin_
        );
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::temperatureCelcius::temperatureCelcius
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "temperatureCelcius", "T"),
    zero_celsius_in_kelvin_
    (
        "zero_celsius_in_kelvin",
        dimTemperature,
        scalar(273.15)
    )
{
    if(!dict.found("result"))
    {
        resultName_ = "Temperature";
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::temperatureCelcius::~temperatureCelcius()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::temperatureCelcius::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    return true;
}

// ************************************************************************* //
