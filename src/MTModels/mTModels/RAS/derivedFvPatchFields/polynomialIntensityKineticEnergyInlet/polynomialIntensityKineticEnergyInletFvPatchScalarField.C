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

#include "polynomialIntensityKineticEnergyInletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polynomialIntensityKineticEnergyInletFvPatchScalarField::
polynomialIntensityKineticEnergyInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    z_(pTraits<vector>::zero),
    intensity_(),
    UName_("U"),
    zGround_(0.0)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


Foam::polynomialIntensityKineticEnergyInletFvPatchScalarField::
polynomialIntensityKineticEnergyInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
    z_(dict.lookup("z")),
    intensity_(Function1<scalar>::New("intensity", dict)),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    zGround_("zGround", dict, p.size())
{
    this->phiName_ = dict.lookupOrDefault<word>("phi","phi");

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


Foam::polynomialIntensityKineticEnergyInletFvPatchScalarField::
polynomialIntensityKineticEnergyInletFvPatchScalarField
(
    const polynomialIntensityKineticEnergyInletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(ptf, p, iF, mapper),
    z_(ptf.z_),
    intensity_(ptf.intensity_, false),
    UName_(ptf.UName_),
    zGround_(mapper(ptf.zGround_))
{}


Foam::polynomialIntensityKineticEnergyInletFvPatchScalarField::
polynomialIntensityKineticEnergyInletFvPatchScalarField
(
    const polynomialIntensityKineticEnergyInletFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(ptf, iF),
    z_(ptf.z_),
    intensity_(ptf.intensity_, false),
    UName_(ptf.UName_),
    zGround_(ptf.zGround_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polynomialIntensityKineticEnergyInletFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    inletOutletFvPatchScalarField::autoMap(m);

    m(zGround_, zGround_);
}


void Foam::polynomialIntensityKineticEnergyInletFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    inletOutletFvPatchScalarField::rmap(ptf, addr);

    const polynomialIntensityKineticEnergyInletFvPatchScalarField& tiptf =
        refCast<const polynomialIntensityKineticEnergyInletFvPatchScalarField>
        (ptf);

    zGround_.rmap(tiptf.zGround_, addr);
}


void Foam::polynomialIntensityKineticEnergyInletFvPatchScalarField::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const vectorField& c = patch().Cf();
    const scalarField coord(c & z_);
    scalarField intensityField(coord.size());

    const fvPatchVectorField& Up =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    const fvsPatchScalarField& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(this->phiName_);

    forAll(coord, i)
    {
      intensityField[i] = intensity_->value(coord[i]-zGround_[i]);
    }

    this->refValue() = 1.5*sqr(intensityField)*magSqr(Up);
    this->valueFraction() = 1.0 - pos0(phip);

    inletOutletFvPatchScalarField::updateCoeffs();
}

void Foam::polynomialIntensityKineticEnergyInletFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "z", z_);
    writeEntry(os, intensity_());
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", this->phiName_);
    writeEntry(os, "zGround", zGround_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        polynomialIntensityKineticEnergyInletFvPatchScalarField
    );
}

// ************************************************************************* //
