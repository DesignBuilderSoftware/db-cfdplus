/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2015-2018 FSD blueCAPE Lda  http://www.bluecape.com.pt
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

#include "uniformJumpPressureMasterFvPatchScalarField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::uniformJumpPressureMasterFvPatchScalarField::applyJump
(
    const scalar & localAverage,
    const scalar & neighbourAverage,
    const scalar & volumetricFlowRate
)
{
    // a simpler way of doing this would be nice
    const scalar specifiedJump = userSpecifiedJump_->value(-volumetricFlowRate);

    this->operator==
    (
        relaxation_*(neighbourAverage + specifiedJump)
      + (1-relaxation_)*localAverage
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniformJumpPressureMasterFvPatchScalarField::
uniformJumpPressureMasterFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    uniformJumpPressureFvPatchField<scalar>(p, iF)
{}


Foam::uniformJumpPressureMasterFvPatchScalarField::
uniformJumpPressureMasterFvPatchScalarField
(
    const uniformJumpPressureMasterFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    uniformJumpPressureFvPatchField<scalar>(ptf, p, iF, mapper)
{}


Foam::uniformJumpPressureMasterFvPatchScalarField::
uniformJumpPressureMasterFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    uniformJumpPressureFvPatchField<scalar>(p, iF, dict)
{}


Foam::uniformJumpPressureMasterFvPatchScalarField::
uniformJumpPressureMasterFvPatchScalarField
(
    const uniformJumpPressureMasterFvPatchScalarField& ptf
)
:
    uniformJumpPressureFvPatchField<scalar>(ptf)
{}


Foam::uniformJumpPressureMasterFvPatchScalarField::
uniformJumpPressureMasterFvPatchScalarField
(
    const uniformJumpPressureMasterFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    uniformJumpPressureFvPatchField<scalar>(ptf, iF)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{
   makePatchTypeField
   (
       fvPatchScalarField,
       uniformJumpPressureMasterFvPatchScalarField
   );
}

// ************************************************************************* //
