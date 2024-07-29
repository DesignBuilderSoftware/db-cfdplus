/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Derived from OpenFOAM's scalarTransport function object.
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

#include "estimateTemperatureOnHeatFluxPatches.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(estimateTemperatureOnHeatFluxPatches, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        estimateTemperatureOnHeatFluxPatches,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::estimateTemperatureOnHeatFluxPatches::correct()
{
    if
    (
        foundObject<volScalarField>(TinputName_)
     && foundObject<volScalarField>(ToutputName_)
    )
    {
        volScalarField& Tin = lookupObjectRef<volScalarField>(TinputName_);
        volScalarField& Tout = lookupObjectRef<volScalarField>(ToutputName_);

        Info<< "MinMax " << TinputName_ << " before heat flux corrections: "
            << min(Tin).value() << ","
            << max(Tin).value()
            << endl;

        Tout == Tin;
        Tout.correctBoundaryConditions();

        Info<< "MinMax " << ToutputName_ << " after heat flux corrections: "
            << min(Tout).value() << ","
            << max(Tout).value()
            << endl;

        Tout = min(Tout, Tmax_);
        Tout = max(Tout, Tmin_);

        Info<< "MinMax " << ToutputName_ << " after imposing limits: "
            << min(Tout).value() << ","
            << max(Tout).value()
            << endl;

    }
    else
    {
        FatalErrorInFunction
            << "One or both fields not found."
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::estimateTemperatureOnHeatFluxPatches::
estimateTemperatureOnHeatFluxPatches
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    TinputName_(""),
    ToutputName_(""),
    Tmin_("Tmin", dimTemperature, 0),
    Tmax_("Tmax", dimTemperature, great)
{
    read(dict);

    if(onConstruction_)
    {
        correct();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::estimateTemperatureOnHeatFluxPatches::
~estimateTemperatureOnHeatFluxPatches()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::estimateTemperatureOnHeatFluxPatches::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);

    onConstruction_ = dict.lookup<word>("onConstruction");
    neverExecute_ = dict.lookup<word>("neverExecute");

    TinputName_ = dict.lookupOrDefault<word>("Tinput", "T");
    ToutputName_ = dict.lookup<word>("Toutput");

    Tmin_ = dimensionedScalar::lookupOrDefault
    (
        "Tmin",
        dict,
        dimTemperature,
        0.0
    );

    Tmax_ = dimensionedScalar::lookupOrDefault
    (
        "Tmax",
        dict,
        dimTemperature,
        great
    );

    if(debug)
    {
        Info << "Tmin: " << Tmin_ << ", Tmax: " << Tmax_ << endl;
    }

    return true;
}


bool Foam::functionObjects::estimateTemperatureOnHeatFluxPatches::execute()
{
    if(!neverExecute_)
    {
        correct();
    }

    return true;
}


bool Foam::functionObjects::estimateTemperatureOnHeatFluxPatches::write()
{
    return true;
}


// ************************************************************************* //
