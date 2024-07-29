/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2015-2020 FSD blueCAPE Lda  http://www.bluecape.com.pt
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

#include "fanRedirectedOutletVelocityFvPatchField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fanRedirectedOutletVelocityFvPatchField::fanRedirectedOutletVelocityFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    remotePatchName_(),
    UName_("U"),
    phiName_("phi"),
    rhoName_("rho"),
    autoAdjustable_(false),
    inletDir_(p.size(), Zero),
    direction_(p.size(), Zero)
{}


fanRedirectedOutletVelocityFvPatchField::fanRedirectedOutletVelocityFvPatchField
(
    const fanRedirectedOutletVelocityFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    remotePatchName_(ptf.remotePatchName_),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    autoAdjustable_(ptf.autoAdjustable_),
    inletDir_(mapper(ptf.inletDir_)),
    direction_(mapper(ptf.direction_))
{}


fanRedirectedOutletVelocityFvPatchField::fanRedirectedOutletVelocityFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false),
    remotePatchName_(dict.lookup("remotePatchName")),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "none")),
    autoAdjustable_(dict.lookupOrDefault<Switch>("autoAdjustable", false)),
    inletDir_(p.size(), vector::zero),
    direction_(p.size(), vector::zero)
{
    direction_ = vectorField("direction", dict, p.size());

    direction_ = direction_/mag(direction_);

    inletDir_ = direction_ / ((-this->patch().nf()) & direction_);

    if (dict.found("value"))
    {
        this->operator==(vectorField("value", dict, p.size()));
    }
    else
    {
        this->operator==(vector::zero);
        evaluate(Pstream::commsTypes::blocking);
    }
}


fanRedirectedOutletVelocityFvPatchField::fanRedirectedOutletVelocityFvPatchField
(
    const fanRedirectedOutletVelocityFvPatchField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    remotePatchName_(ptf.remotePatchName_),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    autoAdjustable_(ptf.autoAdjustable_),
    inletDir_(ptf.inletDir_),
    direction_(ptf.direction_)
{}


fanRedirectedOutletVelocityFvPatchField::fanRedirectedOutletVelocityFvPatchField
(
    const fanRedirectedOutletVelocityFvPatchField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    remotePatchName_(ptf.remotePatchName_),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    autoAdjustable_(ptf.autoAdjustable_),
    inletDir_(ptf.inletDir_),
    direction_(ptf.direction_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fanRedirectedOutletVelocityFvPatchField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField::autoMap(m);

    m(inletDir_,inletDir_);
    m(direction_,direction_);
}


void fanRedirectedOutletVelocityFvPatchField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField::rmap(ptf, addr);

    const fanRedirectedOutletVelocityFvPatchField& tiptf =
        refCast<const fanRedirectedOutletVelocityFvPatchField>
        (ptf);

    inletDir_.rmap(tiptf.inletDir_, addr);
    direction_.rmap(tiptf.direction_, addr);
}


void fanRedirectedOutletVelocityFvPatchField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

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
    scalar sumRemotePatchPhi = 0.0;

    if (phiName_ != "none")
    {
        const surfaceScalarField& phi =
            this->db().objectRegistry::lookupObject<surfaceScalarField>(phiName_);
        const scalarField& remotePatchPhi = phi.boundaryField()[remotePatchID];

        if (rhoName_ != "none")
        {
            const volScalarField& rho =
                this->db().objectRegistry::lookupObject<volScalarField>(rhoName_);
            const scalarField& remotePatchRho =
                rho.boundaryField()[remotePatchID];

            sumRemotePatchPhi = gSum(remotePatchPhi/remotePatchRho);
        }
        else
        {
            sumRemotePatchPhi = gSum(remotePatchPhi);
        }
    }
    else if (UName_ != "none")
    {
        //phi is not always initialized properly, so we have to enforce the
        //direct use of the U field instead
        volVectorField& U =
            this->db().objectRegistry::lookupObjectRef<volVectorField>(UName_);

        if (U.boundaryFieldRef()(remotePatchID))
        {
            fvPatchVectorField& remotePatchU =
                U.boundaryFieldRef()[remotePatchID];

            // Ensure that the corresponding outlet velocity patch field is
            // up-to-date
            remotePatchU.evaluate(Pstream::commsTypes::blocking);

            sumRemotePatchPhi = gSum(remotePatch.Sf() & remotePatchU);
        }
    }
    else
    {
        FatalErrorInFunction
            << "The phi and U entries cannot be both none."
            << nl << exit(FatalError);
    }

    const scalar totArea = gSum(p.magSf());
    const scalar avgURemote = sumRemotePatchPhi/totArea;

    this->operator==(avgURemote*inletDir_);

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void fanRedirectedOutletVelocityFvPatchField::write(Ostream& os) const
{
    // Note: do not write value
    fvPatchField<vector>::write(os);

    writeEntry(os, "remotePatchName", remotePatchName_);

    writeEntryIfDifferent<Switch>(os, "autoAdjustable", "false", autoAdjustable_);

    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "none", rhoName_);

    writeEntry(os, "value", *this);
    writeEntry(os, "direction", direction_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
   fvPatchVectorField,
   fanRedirectedOutletVelocityFvPatchField
);

} // End namespace Foam

// ************************************************************************* //
