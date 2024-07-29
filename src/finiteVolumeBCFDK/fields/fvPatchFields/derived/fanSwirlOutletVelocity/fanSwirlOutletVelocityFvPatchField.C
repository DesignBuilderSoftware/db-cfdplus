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

#include "fanSwirlOutletVelocityFvPatchField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fanSwirlOutletVelocityFvPatchField::fanSwirlOutletVelocityFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    remotePatchName_(),
    phiName_("phi"),
    rhoName_("rho"),
    rpm_(0)
{}


fanSwirlOutletVelocityFvPatchField::fanSwirlOutletVelocityFvPatchField
(
    const fanSwirlOutletVelocityFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    remotePatchName_(ptf.remotePatchName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    rpm_(ptf.rpm_)
{}


fanSwirlOutletVelocityFvPatchField::fanSwirlOutletVelocityFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    remotePatchName_(dict.lookup("remotePatchName")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    rpm_(dict.lookup<scalar>("rpm"))
{
    //Cannot evaluate on this construction, so might as well load the value
    fvPatchField<vector>::operator==(vectorField("value", dict, p.size()));
}


fanSwirlOutletVelocityFvPatchField::fanSwirlOutletVelocityFvPatchField
(
    const fanSwirlOutletVelocityFvPatchField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    rpm_(ptf.rpm_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fanSwirlOutletVelocityFvPatchField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const fvPatch& p = this->patch();

    tmp<vectorField> tnPF = p.nf();
    const vectorField& nPF = tnPF();

    label remotePatchID =
        p.patch().boundaryMesh().findPatchID(remotePatchName_);

    if (remotePatchID < 0)
    {
        FatalErrorInFunction
            << "Unable to find remote patch " << remotePatchName_
            << abort(FatalError);
    }

    const fvPatch& remotePatch = p.boundaryMesh()[remotePatchID];

    const surfaceScalarField& phi =
        this->db().objectRegistry::lookupObject<surfaceScalarField>(phiName_);

    const scalarField& remotePatchPhi = phi.boundaryField()[remotePatchID];
    scalar sumRemotePatchPhi = gSum(remotePatchPhi);

    const scalar totArea = gSum(remotePatch.magSf());
    const scalar avgURemote = sumRemotePatchPhi/totArea;

    const vector avgCenter = gSum(p.Cf() * p.magSf())/totArea;
    const vector avgNormal = -gSum(p.Sf())/totArea;

    // Update angular velocity - convert [rpm] to [rad/s]
    vectorField tangentialVelocity
    (
        (rpm_*constant::mathematical::pi/30.0)
      * (p.Cf() - avgCenter) ^ avgNormal
    );

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        // volumetric flow-rate
        this->operator==(tangentialVelocity + (-nPF)*avgURemote);
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        // mass flow-rate
        operator==(tangentialVelocity + (-nPF)*avgURemote/rhop);
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of " << phiName_ << " are incorrect" << nl
            << "    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << nl << exit(FatalError);
    }

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void fanSwirlOutletVelocityFvPatchField::write(Ostream& os) const
{
    // Note: do not write value
    fvPatchField<vector>::write(os);

    writeEntry(os, "remotePatchName", remotePatchName_);

    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);

    writeEntry(os, "rpm", rpm_);

    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
   fvPatchVectorField,
   fanSwirlOutletVelocityFvPatchField
);

} // End namespace Foam

// ************************************************************************* //
