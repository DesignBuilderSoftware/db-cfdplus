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

#include "uniformJumpPressureFvPatchField.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class Type>
void Foam::uniformJumpPressureFvPatchField<Type>::applyJump
(
    const Type & localAverage,
    const Type & neighbourAverage,
    const scalar & volumetricFlowRate
)
{
    // a simpler way of doing this would be nice
    const Type specifiedJump = userSpecifiedJump_->value(volumetricFlowRate);

    this->operator==
    (
        relaxation_*(neighbourAverage - specifiedJump)
      + (1-relaxation_)*localAverage
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
uniformJumpPressureFvPatchField<Type>::uniformJumpPressureFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    remotePatchName_(),
    phiName_("phi"),
    rhoName_("rho"),
    userSpecifiedJump_(),
    relaxation_(0.5)
{}


template<class Type>
uniformJumpPressureFvPatchField<Type>::uniformJumpPressureFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& fld
)
:
    fixedValueFvPatchField<Type>(p, iF, fld),
    remotePatchName_(),
    phiName_("phi"),
    rhoName_("rho"),
    userSpecifiedJump_(),
    relaxation_(0.5)
{}


template<class Type>
uniformJumpPressureFvPatchField<Type>::uniformJumpPressureFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF),
    remotePatchName_(dict.lookup("remotePatchName")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    userSpecifiedJump_(Function1<Type>::New("specifiedJump", dict)),
    relaxation_(dict.lookup<scalar>("relaxation"))
{
    //Cannot evaluate on this construction, so might as well load the value
    fvPatchField<Type>::operator==(Field<Type>("value", dict, p.size()));
}


template<class Type>
uniformJumpPressureFvPatchField<Type>::uniformJumpPressureFvPatchField
(
    const uniformJumpPressureFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    remotePatchName_(ptf.remotePatchName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    userSpecifiedJump_(ptf.userSpecifiedJump_, false),
    relaxation_(ptf.relaxation_)
{}


template<class Type>
uniformJumpPressureFvPatchField<Type>::uniformJumpPressureFvPatchField
(
    const uniformJumpPressureFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    remotePatchName_(ptf.remotePatchName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    userSpecifiedJump_(ptf.userSpecifiedJump_, false),
    relaxation_(ptf.relaxation_)
{}


template<class Type>
uniformJumpPressureFvPatchField<Type>::uniformJumpPressureFvPatchField
(
    const uniformJumpPressureFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    remotePatchName_(ptf.remotePatchName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    userSpecifiedJump_(ptf.userSpecifiedJump_, false),
    relaxation_(ptf.relaxation_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void uniformJumpPressureFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const GeometricField<Type, fvPatchField, volMesh>& f
    (
        dynamic_cast<const GeometricField<Type, fvPatchField, volMesh>&>
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

    const fvPatchField<Type>& remotePatchField =
        f.boundaryField()[remotePatchID];

    const surfaceScalarField& phi =
        this->db().objectRegistry::template lookupObject<surfaceScalarField>
            (phiName_);

    scalarField remotePatchPhi(phi.boundaryField()[remotePatchID]);

    scalarField localPatchPhi
    (
        this->patch().template lookupPatchField<surfaceScalarField, scalar>
            (phiName_)
    );

    if (rhoName_ != "none")
    {
        const volScalarField& rho =
            this->db().objectRegistry::template lookupObject<volScalarField>
                (rhoName_);

        remotePatchPhi /= rho.boundaryField()[remotePatchID];

        localPatchPhi /=
            this->patch().template lookupPatchField<volScalarField, scalar>
                (rhoName_);
    }

    scalar sumRemotePatchPhi = gSum(remotePatchPhi);
    scalar sumLocalPatchPhi = gSum(localPatchPhi);

    if (mag(sumRemotePatchPhi) > SMALL)
    {
        Type averageRemoteField =
            gSum(remotePatchPhi*remotePatchField)
          / sumRemotePatchPhi;

        Type averageLocalField =
            gSum(*this * localPatchPhi)
           /sumLocalPatchPhi;

        this->applyJump(
            averageLocalField,
            averageRemoteField,
            sumLocalPatchPhi
            );
    }
    else
    {
        Type averageRemoteField =
            gSum(remotePatch.magSf()*remotePatchField)
           /gSum(remotePatch.magSf());

        Type averageLocalField =
            gSum(*this * this->patch().magSf())
           /gSum(this->patch().magSf());

        this->applyJump(
            averageLocalField,
            averageRemoteField,
            sumLocalPatchPhi
            );
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void uniformJumpPressureFvPatchField<Type>::write(Ostream& os) const
{
    // Note: do not write value
    fvPatchField<Type>::write(os);

    writeEntry(os, "remotePatchName", remotePatchName_);

    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);

    writeEntry(os, userSpecifiedJump_());

    writeEntry(os, "relaxation", relaxation_);

    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
