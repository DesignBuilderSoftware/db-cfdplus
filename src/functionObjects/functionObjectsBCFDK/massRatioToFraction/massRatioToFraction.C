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

#include "massRatioToFraction.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(massRatioToFraction, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        massRatioToFraction,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::massRatioToFraction::calc()
{
    if
    (
        foundObject<volScalarField>(fieldName_)
    )
    {
        const volScalarField& field = lookupObject<volScalarField>(fieldName_);

        return store
        (
            resultName_,
            field/(dimensionedScalar("1",field.dimensions(),1) + field)
        );
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::massRatioToFraction::massRatioToFraction
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "massRatioToFraction", "Scalar")
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::massRatioToFraction::~massRatioToFraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::massRatioToFraction::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    return true;
}

// ************************************************************************* //
