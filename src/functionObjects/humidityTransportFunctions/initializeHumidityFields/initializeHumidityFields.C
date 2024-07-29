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

#include "initializeHumidityFields.H"
#include "surfaceFields.H"
#include "dimensionedScalar.H"
#include "addToRunTimeSelectionTable.H"

#include "Common/humidityCalculations.H"
#include "mixedFvPatchField.H"
#include "fixedProfileFvPatchField.H"
#include "zeroGradientFvPatchField.H"
#include "scalarFluxFickFirstLawFvPatchScalarFieldBCFDK.H"
#include "scalarMassFlowRateInletOutletToRatioFvPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(initializeHumidityFields, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        initializeHumidityFields,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::functionObjects::initializeHumidityFields::densityField
(
    const Foam::surfaceScalarField& phi,
    const Foam::dictionary& transportProperties
) const
{
    if (phi.dimensions() == dimMass/dimTime)
    {
        Info<< "Retrieving field " << rhoName_ << endl;
        return mesh_.lookupObject<volScalarField>(rhoName_);
    }
    else if (phi.dimensions() == dimVolume/dimTime)
    {
        const dimensionedScalar rhoRef("rho", dimDensity, transportProperties);

        Info<< "Creating field " << rhoName_ << endl;
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                rhoRef
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "Incompatible dimensions for phi: " << phi.dimensions()
            << nl << "Dimensions should be " << dimMass/dimTime << " or "
            << dimVolume/dimTime << exit(FatalError);
        return tmp<volScalarField>(nullptr);
    }
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::initializeHumidityFields::temperatureField
(
    const Foam::surfaceScalarField& phi,
    const Foam::dictionary& transportProperties
) const
{
    if (phi.dimensions() == dimMass/dimTime)
    {
        Info<< "Retrieving field " << tempName_ << endl;
        return mesh_.lookupObject<volScalarField>(tempName_);
    }
    else if (phi.dimensions() == dimVolume/dimTime)
    {
        dimensionedScalar TRef
        ("TRef", dimTemperature, transportProperties);

        Info<< "Creating field " << tempName_ << endl;
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    tempName_,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                TRef
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "Incompatible dimensions for phi: " << phi.dimensions()
            << nl << "Dimensions should be " << dimMass/dimTime << " or "
            << dimVolume/dimTime << exit(FatalError);
        return tmp<volScalarField>(nullptr);
    }
}


bool Foam::functionObjects::initializeHumidityFields::calc()
{
    bool finalStatus = false;

    if (foundObject<volScalarField>(fieldName_))
    {
        // Create/load the necessary fields and parameters --------------------

        const volScalarField& phi_H2O =
            lookupObject<volScalarField>(fieldName_);

        const surfaceScalarField& phi =
            mesh_.lookupObject<surfaceScalarField>(phiName_);

        const volScalarField& pressure =
            mesh_.lookupObject<volScalarField>(pressureName_);

        dictionary& transportProperties =
            obr_.lookupObjectRef<dictionary>("transportProperties");


        const dimensionedScalar RAir
        ("RAir", dimGasConstant, transportProperties);

        const dimensionedScalar RWaterVapour
        ("RWaterVapour", dimGasConstant, transportProperties);

        const dimensionedScalar PRef
        ("PRef", dimPressure, transportProperties);


        tmp<volScalarField> tT(temperatureField(phi, transportProperties));
        const volScalarField& temperature = tT();

        tmp<volScalarField> tRho(densityField(phi, transportProperties));
        const volScalarField& rho = tRho();

        Info<< "Creating field airPressure" << endl;
        finalStatus = store
        (
            tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "airPressure",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    phi.dimensions() == dimMass/dimTime
                  ?
                    pressure + PRef
                  : // phi.dimensions() == dimVolume/dimTime
                    pressure*rho + PRef
                )
            )
        );

        const volScalarField& airPressure =
            mesh_.lookupObject<volScalarField>("airPressure");

        // Create/load the main fields ----------------------------------------

        Info<< "Creating field absolute humidity" << endl;

        tmp<volScalarField> tw_H2O
        (
            new volScalarField
            (
                IOobject
                (
                    "w_H2O",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                phi_H2O
            )
        );
        volScalarField& w_H2O = tw_H2O.ref();
        w_H2O.dimensions().reset(dimless);

        const volScalarField::Boundary& patchPhi_H2O = phi_H2O.boundaryField();
        volScalarField::Boundary& patch_w_H2O = w_H2O.boundaryFieldRef();

        forAll(patchPhi_H2O, patchi)
        {
            if(isA<fixedProfileFvPatchField<scalar>>(patchPhi_H2O[patchi]))
            {
                patch_w_H2O.set
                (
                    patchi,
                    fvPatchScalarField::New
                    (
                        "fixedValue",
                        patchPhi_H2O[patchi].patch(),
                        w_H2O
                    )
                );
            }
        }


        Info<< "Creating field humidity mass fraction" << endl;
        tmp<volScalarField> tx_H2O
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
                w_H2O
            )
        );
        volScalarField& x_H2O = tx_H2O.ref();
        x_H2O.dimensions().reset(dimless);

        const dimensionedScalar RAirOverRH2O("RAirOverRH2O", RAir/RWaterVapour);
        transportProperties.set("RAirOverRH2O", RAirOverRH2O);

        const volScalarField dimlessTemperature
        (
            temperature
          * dimensionedScalar("invT", dimless/dimTemperature, 1.0)
        );

        Info<< "Creating field saturationPressure" << endl;
        tmp<volScalarField> tsaturationPressure
        (
            new volScalarField
            (
                IOobject
                (
                    "saturationPressure",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                calcSaturationPressure(dimlessTemperature)
              * dimensionedScalar("pressureDim", dimPressure, 1.0)
            )
        );
        const volScalarField& saturationPressure = tsaturationPressure();

        Info<< "Calculating humidity mass ratio" << endl;
        calcHumidityFraction
        (
            w_H2O,
            RAirOverRH2O,
            phi_H2O,
            saturationPressure,
            airPressure
        );

        Info<< "Calculating humidity mass fraction" << endl;
        calcMassFraction(x_H2O, w_H2O);

        Info<< "Calculating humidity fractions for boundaries" << endl;
        {
            const volScalarField::Boundary& patchSaturationPressure =
                saturationPressure.boundaryField();
            const volScalarField::Boundary& patchAirPressure =
                airPressure.boundaryField();
            volScalarField::Boundary& patch_x_H2O = x_H2O.boundaryFieldRef();

            forAll(patchPhi_H2O, patchi)
            {

                if
                (
                    !patchPhi_H2O[patchi].coupled() &&
                    !isA<zeroGradientFvPatchField<scalar>>
                    (
                        patchPhi_H2O[patchi]
                    )
                 &&
                    !isA<scalarFluxFickFirstLawFvPatchScalarFieldBCFDK>
                    (
                        patchPhi_H2O[patchi]
                    )
                 &&
                    !isA<scalarMassFlowRateInletOutletToRatioFvPatchField>
                    (
                        patchPhi_H2O[patchi]
                    )
                )
                {
                    calcHumidityFraction
                    (
                        patch_w_H2O[patchi],
                        RAirOverRH2O.value(),
                        patchPhi_H2O[patchi],
                        patchSaturationPressure[patchi],
                        patchAirPressure[patchi]
                    );

                    calcMassFraction(patch_x_H2O[patchi], patch_w_H2O[patchi]);
                }

                if(isA<mixedFvPatchField<scalar>>(patchPhi_H2O[patchi]))
                {
                    mixedFvPatchField<scalar>& patch_omega_H2O =
                        refCast<mixedFvPatchField<scalar>>(patch_w_H2O[patchi]);

                    mixedFvPatchField<scalar>& patch_mFrac_H2O =
                        refCast<mixedFvPatchField<scalar>>(patch_x_H2O[patchi]);

                    patch_omega_H2O.refValue() = patch_omega_H2O;
                    patch_mFrac_H2O.refValue() = patch_mFrac_H2O;
                }
            }

        }

        Info<< "Calculating maximum absolute humidity" << endl;
        {
            const volScalarField phi_H2O_max
            (
                IOobject
                (
                  "phi_H2O_max",
                  mesh_.time().timeName(),
                  mesh_,
                  IOobject::NO_READ,
                  IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("phi_H2O_max",dimless,1.0)
            );

            volScalarField max_w_H2O
            (
              IOobject
              (
                  "max_w_H2O",
                  mesh_.time().timeName(),
                  mesh_,
                  IOobject::NO_READ,
                  IOobject::NO_WRITE
              ),
              phi_H2O_max
            );

            calcHumidityFraction
            (
                max_w_H2O,
                RAirOverRH2O,
                phi_H2O_max,
                saturationPressure,
                airPressure
            );

            tmp<volScalarField> twr_H2O
            (
                new volScalarField
                (
                    IOobject
                    (
                        "maxAbsoluteHumidity",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    max_w_H2O * rho
                )
            );

            Info<< "Calculating maximum humidity mass fraction" << endl;
            tmp<volScalarField> tmax_x_H2O
            (
                new volScalarField
                (
                    IOobject
                    (
                        "max_"+resultName_,
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    x_H2O
                )
            );
            volScalarField& max_x_H2O = tmax_x_H2O.ref();

            calcMassFraction(max_x_H2O, max_w_H2O);

            finalStatus = store(twr_H2O) && finalStatus;
            finalStatus = store(tmax_x_H2O) && finalStatus;
        }

        finalStatus = store(tw_H2O) && finalStatus;
        finalStatus = store(tx_H2O) && finalStatus;
        finalStatus = store(tsaturationPressure) && finalStatus;
        finalStatus = store(tRho) && finalStatus;

    }

    return finalStatus;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::initializeHumidityFields::initializeHumidityFields
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "initializeHumidityFields", "phi_H2O")
{
    resultName_ = "x_H2O";
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::initializeHumidityFields::~initializeHumidityFields()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::initializeHumidityFields::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    debug = dict.lookupOrDefault<label>("debug", debug);

    phiName_ = dict.lookupOrDefault<word>("phi", "phi");
    rhoName_ = dict.lookupOrDefault<word>("rho", "rho");
    tempName_ = dict.lookupOrDefault<word>("temperature", "T");
    pressureName_ = dict.lookupOrDefault<word>("pressure", "p");

    return true;
}

bool Foam::functionObjects::initializeHumidityFields::write()
{
    bool status = writeObject(resultName_);

    if(debug)
    {
        writeObject("w_H2O");
        writeObject("max_"+resultName_);
    }

    status = writeObject("maxAbsoluteHumidity") && status;

    return status;
}

// ************************************************************************* //
