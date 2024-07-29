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

#include "limitedTurbulentMixingLengthFrequencyInletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "momentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

limitedTurbulentMixingLengthFrequencyInletFvPatchScalarField::
limitedTurbulentMixingLengthFrequencyInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    mixingLength_(0.0),
    kName_("undefined-k")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}

limitedTurbulentMixingLengthFrequencyInletFvPatchScalarField::
limitedTurbulentMixingLengthFrequencyInletFvPatchScalarField
(
    const limitedTurbulentMixingLengthFrequencyInletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(ptf, p, iF, mapper),
    mixingLength_(ptf.mixingLength_),
    kName_(ptf.kName_)
{}

limitedTurbulentMixingLengthFrequencyInletFvPatchScalarField::
limitedTurbulentMixingLengthFrequencyInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
    mixingLength_(dict.lookup<scalar>("mixingLength")),
    kName_(dict.lookupOrDefault<word>("k", "k"))
{
    this->phiName_ = dict.lookupOrDefault<word>("phi", "phi");

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}

limitedTurbulentMixingLengthFrequencyInletFvPatchScalarField::
limitedTurbulentMixingLengthFrequencyInletFvPatchScalarField
(
    const limitedTurbulentMixingLengthFrequencyInletFvPatchScalarField& ptf
)
:
    inletOutletFvPatchScalarField(ptf),
    mixingLength_(ptf.mixingLength_),
    kName_(ptf.kName_)
{}

limitedTurbulentMixingLengthFrequencyInletFvPatchScalarField::
limitedTurbulentMixingLengthFrequencyInletFvPatchScalarField
(
    const limitedTurbulentMixingLengthFrequencyInletFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(ptf, iF),
    mixingLength_(ptf.mixingLength_),
    kName_(ptf.kName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void limitedTurbulentMixingLengthFrequencyInletFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Lookup Cmu corresponding to the turbulence model selected
    const momentumTransportModel& turbModel =
        db().lookupObject<momentumTransportModel>
        (
            IOobject::groupName
            (
                momentumTransportModel::typeName,
                internalField().group()
            )
        );

    const scalar Cmu =
        turbModel.coeffDict().lookupOrDefault<scalar>("Cmu", 0.09);

    const scalar Cmu25 = pow(Cmu, 0.25);

    const fvPatchScalarField& kp =
        patch().lookupPatchField<volScalarField, scalar>(kName_);

    const fvsPatchScalarField& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(this->phiName_);

    this->refValue() = max(sqrt(kp)/(Cmu25*mixingLength_), ROOTVSMALL);
    this->valueFraction() = 1.0 - pos0(phip);

    inletOutletFvPatchScalarField::updateCoeffs();
}


void limitedTurbulentMixingLengthFrequencyInletFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "mixingLength", mixingLength_);
    writeEntry(os, "phi", this->phiName_);
    writeEntry(os, "k", kName_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    limitedTurbulentMixingLengthFrequencyInletFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
