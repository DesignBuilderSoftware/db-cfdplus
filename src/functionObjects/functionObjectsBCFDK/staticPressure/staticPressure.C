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

#include "staticPressure.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(staticPressure, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        staticPressure,
        dictionary
    );
}
}


template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::staticPressure::modeEnum,
    2
>::names[] =
{
    "isothermal",
    "buoyant"
};

const Foam::NamedEnum<Foam::functionObjects::staticPressure::modeEnum, 2>
    Foam::functionObjects::staticPressure::modeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::staticPressure::calc()
{
    if
    (
        foundObject<volScalarField>(fieldName_)
    )
    {
        const volScalarField& p = lookupObject<volScalarField>(fieldName_);

        if (mode_ == mdIsothermal)
        {
            return store
            (
                resultName_,
                p*rho_
            );
        }
        else //if(mode_ == mdBuoyant)
        {
            return store
            (
                resultName_,
                p*scalar(1.0) //dummy operation for automatic recasting
            );
        }
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::staticPressure::staticPressure
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "staticPressure", "p"),
    mode_(mdIsothermal)
{
    resultName_ = "Pressure";
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::staticPressure::~staticPressure()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::staticPressure::read(const dictionary& dict)
{
    Info<< type() << " read dictionary:" << endl;

    fieldExpression::read(dict);

    mode_ = modeNames_
    [
        dict.lookupOrDefault<word>("mode", "isothermal")
    ];

    if(mode_ == mdIsothermal)
    {
        const dictionary& transportProperties =
            obr_.lookupObject<dictionary>("transportProperties");

        Info<< "    Reading density rho\n" << endl;

        transportProperties.lookup("rho") >> rho_;
    }

    return true;
}

// ************************************************************************* //
