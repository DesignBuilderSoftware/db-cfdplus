/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2015-2021 FSD blueCAPE Lda  http://www.bluecape.com.pt
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

#include "uniformPressureFunctionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

#define CHANGE_RELAX_FACTOR_BY 0.99

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* NamedEnum
    <
        uniformPressureFunctionFvPatchScalarField::fanFlowDirection,
        2
    >::names[] =
    {
        "in",
        "out"
    };
}

const Foam::NamedEnum
<
    Foam::uniformPressureFunctionFvPatchScalarField::fanFlowDirection,
    2
> Foam::uniformPressureFunctionFvPatchScalarField::fanFlowDirectionNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniformPressureFunctionFvPatchScalarField::
  uniformPressureFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    dynamicPressureFvPatchScalarField(p, iF),
    function_(),
    direction_(ffdOut),
    phiName_("phi"),
    minFlowRate_(-1.0),
    maxFlowRate_(-1.0),
    relaxation_(1.0),
    prevValue_(0),
    previousTimeIndex_(db().time().timeIndex())
{}


Foam::uniformPressureFunctionFvPatchScalarField::
  uniformPressureFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    dynamicPressureFvPatchScalarField(p, iF, dict),
    function_(Function1<scalar>::New("function", dict)),
    direction_(fanFlowDirectionNames_.read(dict.lookup("direction"))),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    minFlowRate_(dict.lookupOrDefault<scalar>("minFlowRate", -1.0)),
    maxFlowRate_(dict.lookupOrDefault<scalar>("maxFlowRate", -1.0)),
    relaxation_(dict.lookupOrDefault<scalar>("relaxation", 1.0)),
    prevValue_(dict.lookupOrDefault<scalar>("prevValue", 0.0)),
    previousTimeIndex_(db().time().timeIndex())
{}


Foam::uniformPressureFunctionFvPatchScalarField::
  uniformPressureFunctionFvPatchScalarField
(
    const uniformPressureFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    dynamicPressureFvPatchScalarField(ptf, p, iF, mapper),
    function_(ptf.function_, false),
    direction_(ptf.direction_),
    phiName_(ptf.phiName_),
    minFlowRate_(ptf.minFlowRate_),
    maxFlowRate_(ptf.maxFlowRate_),
    relaxation_(ptf.relaxation_),
    prevValue_(ptf.prevValue_),
    previousTimeIndex_(ptf.previousTimeIndex_)
{}


Foam::uniformPressureFunctionFvPatchScalarField::
  uniformPressureFunctionFvPatchScalarField
(
    const uniformPressureFunctionFvPatchScalarField& pfopsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    dynamicPressureFvPatchScalarField(pfopsf, iF),
    function_(pfopsf.function_, false),
    direction_(pfopsf.direction_),
    phiName_(pfopsf.phiName_),
    minFlowRate_(pfopsf.minFlowRate_),
    maxFlowRate_(pfopsf.maxFlowRate_),
    relaxation_(pfopsf.relaxation_),
    prevValue_(pfopsf.prevValue_),
    previousTimeIndex_(pfopsf.previousTimeIndex_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uniformPressureFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Retrieve flux field
    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    scalar flowSign = direction_==ffdOut ? +1 : -1;
    scalar pressureSign = direction_==ffdIn ? +1 : -1;

    // Average volumetric flow rate
    //scalar unfiltered = 0;
    scalar volFlowRate = 0;

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        volFlowRate = gSum(max(flowSign*phip,0.0));
        //unfiltered = flowSign*gSum(phip);
    }
    else if (phi.dimensions() == dimVelocity*dimArea*dimDensity)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);
        volFlowRate = gSum(max(flowSign*phip,0.0)/rhop);
        //unfiltered = flowSign*gSum(phip/rhop);
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of phi are not correct"
            << "\n    on patch " << patch().name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath() << nl
            << exit(FatalError);
    }

    if (minFlowRate_ >= 0.0)
    {
        volFlowRate = max(volFlowRate, minFlowRate_);
    }

    if (maxFlowRate_ > 0.0)
    {
        volFlowRate = min(volFlowRate, maxFlowRate_);
    }

    // Pressure drop for this flow rate
    scalar pdFan = function_->value(volFlowRate);

    const scalar targetPd =
        (1-relaxation_) * prevValue_
      + relaxation_ * pressureSign*pdFan;

    operator==(targetPd);

    // Update the old-time quantities
    if (previousTimeIndex_ != db().time().timeIndex())
    {
        previousTimeIndex_ = db().time().timeIndex();
        prevValue_ = targetPd;
    }

    fixedValueFvPatchScalarField::updateCoeffs();

    /*
    const scalar area = gSum(this->patch().magSf());
    const scalar Unp = volFlowRate/area;

    dynamicPressureFvPatchScalarField::updateCoeffs
    (
        p0_,
        scalarField(patch().size(), pressureSign*pdFan)
        - 0.5*Unp*mag(Unp)
    );*/
}


void Foam::uniformPressureFunctionFvPatchScalarField::write(Ostream& os) const
{
    dynamicPressureFvPatchScalarField::write(os);
    writeEntry(os, function_());
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntry(os, "direction", fanFlowDirectionNames_[direction_]);

    if (minFlowRate_ >= 0.0)
    {
        writeEntry(os, "minFlowRate", minFlowRate_);
    }

    if (maxFlowRate_ > 0.0)
    {
        writeEntry(os, "maxFlowRate", maxFlowRate_);
    }

    writeEntryIfDifferent<scalar>(os, "relaxation", 1.0, relaxation_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if defined( WIN32 ) || defined( WIN64 )
#include "PrghPressureFvPatchScalarField.T.H"
#else
#include "PrghPressureFvPatchScalarField.H"
#endif

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        uniformPressureFunctionFvPatchScalarField
    );

    makePrghPatchScalarField
    (
        uniformPressureFunction,
        prghUniformPressureFunction
    )
};


// ************************************************************************* //
