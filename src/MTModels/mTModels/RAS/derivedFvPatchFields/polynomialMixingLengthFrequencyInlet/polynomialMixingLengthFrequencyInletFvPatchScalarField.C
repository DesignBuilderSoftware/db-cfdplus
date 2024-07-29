/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's atmBoundaryLayer* incompressible boundary condition
    classes.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2014-2019 FSD blueCAPE Lda  http://www.bluecape.com.pt
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

#include "polynomialMixingLengthFrequencyInletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "momentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

polynomialMixingLengthFrequencyInletFvPatchScalarField::
polynomialMixingLengthFrequencyInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    z_(pTraits<vector>::zero),
    mixingLength_(),
    kName_("undefined-k"),
    zGround_(0.0)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}

polynomialMixingLengthFrequencyInletFvPatchScalarField::
polynomialMixingLengthFrequencyInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
    z_(dict.lookup("z")),
    mixingLength_(Function1<scalar>::New("mixingLength", dict)),
    kName_(dict.lookupOrDefault<word>("k", "k")),
    zGround_("zGround", dict, p.size())
{
    this->phiName_ = dict.lookupOrDefault<word>("phi", "phi");

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}

polynomialMixingLengthFrequencyInletFvPatchScalarField::
polynomialMixingLengthFrequencyInletFvPatchScalarField
(
    const polynomialMixingLengthFrequencyInletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(ptf, p, iF, mapper),
    z_(ptf.z_),
    mixingLength_(ptf.mixingLength_().clone().ptr()),
    kName_(ptf.kName_),
    zGround_(mapper(ptf.zGround_))
{}

polynomialMixingLengthFrequencyInletFvPatchScalarField::
polynomialMixingLengthFrequencyInletFvPatchScalarField
(
    const polynomialMixingLengthFrequencyInletFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(ptf, iF),
    z_(ptf.z_),
    mixingLength_(ptf.mixingLength_().clone().ptr()),
    kName_(ptf.kName_),
    zGround_(ptf.zGround_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polynomialMixingLengthFrequencyInletFvPatchScalarField::updateCoeffs()
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

    const vectorField& c = patch().Cf();
    const scalarField coord(c & z_);
    scalarField mixingLengthField(coord.size());

    forAll(coord, i)
    {
      mixingLengthField[i] = mixingLength_->value(coord[i]-zGround_[i]);
    }

    this->refValue() = sqrt(kp)/(Cmu25*mixingLengthField);
    this->valueFraction() = 1.0 - pos0(phip);

    inletOutletFvPatchScalarField::updateCoeffs();
}


void polynomialMixingLengthFrequencyInletFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "z", z_);
    writeEntry(os, mixingLength_());
    writeEntry(os, "phi", this->phiName_);
    writeEntry(os, "k", kName_);
    writeEntry(os, "zGround", zGround_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    polynomialMixingLengthFrequencyInletFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
