/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Modified and developed for blueCFD(R)-Kernel
    Derived from blueCFD(R)-Kernel's massAccounting function object.
    Copyright (C) 2020 FSD blueCAPE Lda  http://www.bluecape.com.pt
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

#include "resetBuoyantPressures.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(resetBuoyantPressures, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        resetBuoyantPressures,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::resetBuoyantPressures::rebalance()
{
    volScalarField& p =
        mesh_.lookupObjectRef<volScalarField>("p");

    volScalarField& p_rgh =
        mesh_.lookupObjectRef<volScalarField>("p_rgh");

    const volScalarField& rho =
        mesh_.lookupObject<volScalarField>("rho");

    const volScalarField& gh =
        mesh_.lookupObject<volScalarField>("gh");

    const uniformDimensionedScalarField& pRef =
        lookupObject<uniformDimensionedScalarField>("pRef");

    volScalarField p_rghReRead
    (
        IOobject
        (
            "p_rgh",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_
    );

    //Force re-reading it
    p_rgh = p_rghReRead;

    // Correct the p field accordingly
    p = p_rgh + rho*gh + pRef;

    p.write();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::resetBuoyantPressures::resetBuoyantPressures
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::resetBuoyantPressures::~resetBuoyantPressures()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::resetBuoyantPressures::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    return true;
}


bool Foam::functionObjects::resetBuoyantPressures::execute()
{
    return true;
}


bool Foam::functionObjects::resetBuoyantPressures::write()
{
    rebalance();

    return true;
}


// ************************************************************************* //
