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

#include "finalizeHumidityFields.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "addToRunTimeSelectionTable.H"

#include "Common/humidityCalculations.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(finalizeHumidityFields, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        finalizeHumidityFields,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::finalizeHumidityFields::calc()
{
    bool finalStatus = false;

    if (foundObject<volScalarField>(fieldName_))
    {
        // Create/load the necessary fields and parameters --------------------

        const volScalarField& x_H2O =
            lookupObject<volScalarField>(fieldName_);

        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);

        const dictionary& transportProperties =
            obr_.lookupObject<dictionary>("transportProperties");


        const dimensionedScalar RAir
        ("RAir", dimGasConstant, transportProperties);

        const dimensionedScalar RWaterVapour
        ("RWaterVapour", dimGasConstant, transportProperties);

        const dimensionedScalar RAirOverRH2O
        ("RAirOverRH2O", dimless, transportProperties);


        const volScalarField& airPressure =
            mesh_.lookupObject<volScalarField>("airPressure");

        const volScalarField& saturationPressure =
            mesh_.lookupObject<volScalarField>("saturationPressure");


        // Calculate the relative humidity ------------------------------------

        Info<< "Calculating relative humidity" << endl;

        tmp<volScalarField> tw_H2O_post
        (
            new volScalarField
            (
                IOobject
                (
                  "w_H2O_post",
                  mesh_.time().timeName(),
                  mesh_,
                  IOobject::NO_READ,
                  IOobject::NO_WRITE
                ),
                x_H2O / (scalar(1) - x_H2O)
            )
        );
        const volScalarField& w_H2O_post = tw_H2O_post();

        tmp<volScalarField> trelativeHumidity
        (
            new volScalarField
            (
                IOobject
                (
                  resultName_,
                  mesh_.time().timeName(),
                  mesh_,
                  IOobject::NO_READ,
                  IOobject::NO_WRITE
                ),
                airPressure * w_H2O_post * scalar(100)
              / ((RAirOverRH2O + w_H2O_post) * saturationPressure)
            )
        );
        volScalarField& relativeHumidity = trelativeHumidity.ref();

        if(useMinimumRelativeHumidity_)
        {
            relativeHumidity = max(relativeHumidity, minimumRelativeHumidity_);
        }

        finalStatus = store(trelativeHumidity);

        if (debug)
        {
            finalStatus = store(tw_H2O_post) && finalStatus;
        }


        // Calculate the absolute humidity ------------------------------------

        Info<< "Calculating absolute humidity" << endl;
        finalStatus = store
        (
            "absoluteHumidity",
            w_H2O_post * rho
        ) && finalStatus;
    }

    return finalStatus;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::finalizeHumidityFields::finalizeHumidityFields
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "finalizeHumidityFields", "x_H2O"),
    minimumRelativeHumidity_(-1.0)
{
    resultName_ = "relativeHumidity";
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::finalizeHumidityFields::~finalizeHumidityFields()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::finalizeHumidityFields::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    debug = dict.lookupOrDefault<label>("debug", debug);

    rhoName_ = dict.lookupOrDefault<word>("rho", "rho");

    useMinimumRelativeHumidity_ = dict.found("minimumRelativeHumidity");
    if(useMinimumRelativeHumidity_)
    {
        minimumRelativeHumidity_ =
            dict.lookup<scalar>("minimumRelativeHumidity");
    }

    return true;
}

bool Foam::functionObjects::finalizeHumidityFields::write()
{
    bool status = writeObject(resultName_);

    if(debug)
    {
        writeObject("w_H2O_post");
    }

    status = writeObject("absoluteHumidity") && status;

    return status;
}

// ************************************************************************* //
