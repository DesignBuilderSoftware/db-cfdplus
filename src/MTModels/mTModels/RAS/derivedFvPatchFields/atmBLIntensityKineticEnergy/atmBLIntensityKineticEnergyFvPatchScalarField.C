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

#include "atmBLIntensityKineticEnergyFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::atmBLIntensityKineticEnergyFvPatchScalarField::
atmBLIntensityKineticEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    z_(pTraits<vector>::zero),
    a_(0.22),
    zGround_(0.0),
    UName_("U")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


Foam::atmBLIntensityKineticEnergyFvPatchScalarField::
atmBLIntensityKineticEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
    z_(dict.lookup("z")),
    a_(dict.lookup<scalar>("a")),
    zGround_("zGround", dict, p.size()),
    UName_(dict.lookupOrDefault<word>("U", "U"))
{
    this->phiName_ = dict.lookupOrDefault<word>("phi","phi");

    if (a_ < 0 || a_ > 1)
    {
        FatalErrorInFunction
            << "Inlet exponent should be specified as a fraction 0-1 "
               "of the mean velocity\n"
               "    value given is " << a_ << nl
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


Foam::atmBLIntensityKineticEnergyFvPatchScalarField::
atmBLIntensityKineticEnergyFvPatchScalarField
(
    const atmBLIntensityKineticEnergyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(ptf, p, iF, mapper),
    z_(ptf.z_),
    a_(ptf.a_),
    zGround_(mapper(ptf.zGround_)),
    UName_(ptf.UName_)
{}


Foam::atmBLIntensityKineticEnergyFvPatchScalarField::
atmBLIntensityKineticEnergyFvPatchScalarField
(
    const atmBLIntensityKineticEnergyFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(ptf, iF),
    z_(ptf.z_),
    a_(ptf.a_),
    zGround_(ptf.zGround_),
    UName_(ptf.UName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::atmBLIntensityKineticEnergyFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    inletOutletFvPatchScalarField::autoMap(m);

    m(zGround_, zGround_);
}


void Foam::atmBLIntensityKineticEnergyFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    inletOutletFvPatchScalarField::rmap(ptf, addr);

    const atmBLIntensityKineticEnergyFvPatchScalarField& tiptf =
        refCast<const atmBLIntensityKineticEnergyFvPatchScalarField>
        (ptf);

    zGround_.rmap(tiptf.zGround_, addr);
}


void Foam::atmBLIntensityKineticEnergyFvPatchScalarField::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalar zg = SMALL;
    scalar zb = SMALL;

    if(a_ <= 0.1)
    {
      zg = 250;
      zb = 5;
    }
    else if(0.1 < a_ && a_ <= 0.14)
    {
      zg = 350;
      zb= 5;
    }
    else if(0.14 < a_ && a_ <= 0.22)
    {
      zg = 450;
      zb = 10;
    }
    else
    {
      zg = 650;
      zb = 30;
    }

    const vectorField& c = patch().Cf();
    const scalarField coord(c & z_);
    scalarField intensity(coord.size());

    const fvPatchVectorField& Up =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    const fvsPatchScalarField& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(this->phiName_);


    forAll(coord, i)
    {
      if(coord[i] - zGround_[i] <= zb)
      {
        intensity[i] = 0.1 * pow((zb/zg),-a_-0.05);
      }
      else
      {
        intensity[i] = 0.1 * pow(((coord[i]-zGround_[i])/zg),-a_-0.05);
      }
    }

    this->refValue() = 1.5*sqr(intensity)*magSqr(Up);
    this->valueFraction() = 1.0 - pos0(phip);

    inletOutletFvPatchScalarField::updateCoeffs();
}


void Foam::atmBLIntensityKineticEnergyFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "z", z_);
    writeEntry(os, "a", a_);
    writeEntry(os, "zGround", zGround_);
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
        atmBLIntensityKineticEnergyFvPatchScalarField
    );
}

// ************************************************************************* //
