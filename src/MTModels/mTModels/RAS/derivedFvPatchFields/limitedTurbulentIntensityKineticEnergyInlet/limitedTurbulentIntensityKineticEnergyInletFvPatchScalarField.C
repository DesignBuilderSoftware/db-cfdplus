/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2014-2018 FSD blueCAPE Lda  http://www.bluecape.com.pt
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

#include "limitedTurbulentIntensityKineticEnergyInletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::limitedTurbulentIntensityKineticEnergyInletFvPatchScalarField::
limitedTurbulentIntensityKineticEnergyInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    intensity_(0.0),
    UName_("U")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}

Foam::limitedTurbulentIntensityKineticEnergyInletFvPatchScalarField::
limitedTurbulentIntensityKineticEnergyInletFvPatchScalarField
(
    const limitedTurbulentIntensityKineticEnergyInletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(ptf, p, iF, mapper),
    intensity_(ptf.intensity_),
    UName_(ptf.UName_)
{}

Foam::limitedTurbulentIntensityKineticEnergyInletFvPatchScalarField::
limitedTurbulentIntensityKineticEnergyInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
    intensity_(dict.lookup<scalar>("intensity")),
    UName_(dict.lookupOrDefault<word>("U", "U"))
{
    this->phiName_ = dict.lookupOrDefault<word>("phi", "phi");

    if (intensity_ < 0 || intensity_ > 1)
    {
        FatalErrorInFunction
            << "Turbulence intensity should be specified as a fraction 0-1 "
               "of the mean velocity\n"
               "    value given is " << intensity_ << nl
            << "    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}

Foam::limitedTurbulentIntensityKineticEnergyInletFvPatchScalarField::
limitedTurbulentIntensityKineticEnergyInletFvPatchScalarField
(
    const limitedTurbulentIntensityKineticEnergyInletFvPatchScalarField& ptf
)
:
    inletOutletFvPatchScalarField(ptf),
    intensity_(ptf.intensity_),
    UName_(ptf.UName_)
{}


Foam::limitedTurbulentIntensityKineticEnergyInletFvPatchScalarField::
limitedTurbulentIntensityKineticEnergyInletFvPatchScalarField
(
    const limitedTurbulentIntensityKineticEnergyInletFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(ptf, iF),
    intensity_(ptf.intensity_),
    UName_(ptf.UName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::limitedTurbulentIntensityKineticEnergyInletFvPatchScalarField::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchVectorField& Up =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    const fvsPatchScalarField& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(this->phiName_);

    this->refValue() = max(1.5*sqr(intensity_)*magSqr(Up), ROOTVSMALL);
    this->valueFraction() = 1.0 - pos0(phip);

    inletOutletFvPatchScalarField::updateCoeffs();
}


void Foam::limitedTurbulentIntensityKineticEnergyInletFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "intensity", intensity_);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", this->phiName_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        limitedTurbulentIntensityKineticEnergyInletFvPatchScalarField
    );
}

// ************************************************************************* //
