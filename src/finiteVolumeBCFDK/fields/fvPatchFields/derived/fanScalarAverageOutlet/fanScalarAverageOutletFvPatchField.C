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

#include "fanScalarAverageOutletFvPatchField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fanScalarAverageOutletFvPatchField::fanScalarAverageOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    remotePatchName_(),
    phiName_("phi")
{}


fanScalarAverageOutletFvPatchField::fanScalarAverageOutletFvPatchField
(
    const fanScalarAverageOutletFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    remotePatchName_(ptf.remotePatchName_),
    phiName_(ptf.phiName_)
{}


fanScalarAverageOutletFvPatchField::fanScalarAverageOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    remotePatchName_(dict.lookup("remotePatchName")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi"))
{
    fvPatchField<scalar>::operator==(Field<scalar>("value", dict, p.size()));
}


fanScalarAverageOutletFvPatchField::fanScalarAverageOutletFvPatchField
(
    const fanScalarAverageOutletFvPatchField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    remotePatchName_(ptf.remotePatchName_),
    phiName_(ptf.phiName_)
{}


fanScalarAverageOutletFvPatchField::fanScalarAverageOutletFvPatchField
(
    const fanScalarAverageOutletFvPatchField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    remotePatchName_(ptf.remotePatchName_),
    phiName_(ptf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fanScalarAverageOutletFvPatchField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const GeometricField<scalar, fvPatchField, volMesh>& f
    (
        dynamic_cast<const GeometricField<scalar, fvPatchField, volMesh>&>
        (
            this->internalField()
        )
    );

    const fvPatch& p = this->patch();
    label remotePatchID =
        p.patch().boundaryMesh().findPatchID(remotePatchName_);

    if (remotePatchID < 0)
    {
        FatalErrorInFunction
            << "Unable to find remote patch " << remotePatchName_
            << abort(FatalError);
    }

    const fvPatch& remotePatch = p.boundaryMesh()[remotePatchID];

    const fvPatchField<scalar>& remotePatchField =
        f.boundaryField()[remotePatchID];

    const surfaceScalarField& phi =
        this->db().objectRegistry::lookupObject<surfaceScalarField>(phiName_);

    const scalarField& remotePatchPhi = phi.boundaryField()[remotePatchID];

    const scalar sumRemotePatchPhi = gSum(remotePatchPhi);

    if (sumRemotePatchPhi > SMALL)
    {
        const scalar averageVariableRemote =
            gSum(remotePatchField * remotePatchPhi)
          / sumRemotePatchPhi;

        this->operator==(averageVariableRemote);
    }
    else
    {
        const scalar totArea = gSum(remotePatch.magSf());

        const scalar averageVariableRemote =
            gSum(remotePatchField * remotePatch.magSf())
          / totArea;

        this->operator==(averageVariableRemote);
    }

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void fanScalarAverageOutletFvPatchField::write(Ostream& os) const
{
    // Note: do not write value
    fvPatchField<scalar>::write(os);

    writeEntry(os, "remotePatchName", remotePatchName_);

    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
   fvPatchScalarField,
   fanScalarAverageOutletFvPatchField
);

} // End namespace Foam

// ************************************************************************* //
