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

#include "pressureDirectedSwirlFvPatchField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* NamedEnum
    <
        pressureDirectedSwirlFvPatchField::ventilatorFlowDirection,
        2
    >::names[] =
    {
        "in",
        "out"
    };
}

const Foam::NamedEnum
<
    Foam::pressureDirectedSwirlFvPatchField::ventilatorFlowDirection,
    2
> Foam::pressureDirectedSwirlFvPatchField::ventilatorFlowDirectionNames_;



namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pressureDirectedSwirlFvPatchField::pressureDirectedSwirlFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    phiName_("phi"),
    rhoName_("rho"),
    rpm_(0),
    rpmOrDiffuser_(false),
    minFlowRate_(-1.0),
    maxFlowRate_(-1.0),
    diffuserDir_(p.size(),vector::zero),
    direction_(vfdOut)
{}

pressureDirectedSwirlFvPatchField::pressureDirectedSwirlFvPatchField
(
    const pressureDirectedSwirlFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    rpm_(ptf.rpm_),
    rpmOrDiffuser_(ptf.rpmOrDiffuser_),
    minFlowRate_(ptf.minFlowRate_),
    maxFlowRate_(ptf.maxFlowRate_),
    diffuserDir_(mapper(ptf.diffuserDir_)),
    direction_(ptf.direction_)
{}


pressureDirectedSwirlFvPatchField::pressureDirectedSwirlFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    rpm_(dict.lookupOrDefault<scalar>("rpm", 0.0)),
    rpmOrDiffuser_(dict.found("rpm")),
    minFlowRate_(dict.lookupOrDefault<scalar>("minFlowRate", -1.0)),
    maxFlowRate_(dict.lookupOrDefault<scalar>("maxFlowRate", -1.0)),
    diffuserDir_(p.size(), vector::zero),
    direction_(ventilatorFlowDirectionNames_.read(dict.lookup("direction")))
{
    if (dict.found("diffuser"))
    {
        if (rpmOrDiffuser_)
        {
            FatalErrorInFunction
                << "cannot set both 'rpm' and 'diffuser'" << nl
                << exit(FatalError);
        }

        diffuserDir_ = vectorField("diffuser", dict, p.size());
        diffuserDir_ /= mag(diffuserDir_);
    }

    if (dict.found("value"))
    {
        this->operator==(vectorField("value", dict, p.size()));
    }
    else
    {
        const scalar totArea = gSum(p.magSf());

        if (maxFlowRate_ > 0 && !rpmOrDiffuser_)
        {
          this->operator==(diffuserDir_ * maxFlowRate_/totArea);
        }
        else
        {
          this->operator==(vector::zero);
        }
    }
}


pressureDirectedSwirlFvPatchField::pressureDirectedSwirlFvPatchField
(
    const pressureDirectedSwirlFvPatchField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    rpm_(ptf.rpm_),
    rpmOrDiffuser_(ptf.rpmOrDiffuser_),
    minFlowRate_(ptf.minFlowRate_),
    maxFlowRate_(ptf.maxFlowRate_),
    diffuserDir_(ptf.diffuserDir_),
    direction_(ptf.direction_)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pressureDirectedSwirlFvPatchField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField::autoMap(m);

    m(diffuserDir_,diffuserDir_);
}


void pressureDirectedSwirlFvPatchField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField::rmap(ptf, addr);

    const pressureDirectedSwirlFvPatchField& tiptf =
        refCast<const pressureDirectedSwirlFvPatchField>
        (ptf);

    diffuserDir_.rmap(tiptf.diffuserDir_, addr);
}


void pressureDirectedSwirlFvPatchField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const fvPatch& p = this->patch();

    vectorField nPF(p.nf());

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    scalar sumPatchPhi = 0.0;

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        // volumetric flow-rate
        sumPatchPhi = gSum(phip);
    }
    else if (phi.dimensions() == dimVelocity*dimArea*dimDensity)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        // mass flow-rate
        sumPatchPhi = gSum(phip/rhop);
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

    int dir = 2*direction_ - 1;

    const scalar totArea = gSum(p.magSf());
    const scalar avgURemote = dir*sumPatchPhi/totArea;

    vectorField tangentialVelocity(p.size(), vector::zero);
    scalarField normalVelocity(p.size(), 0);

    if(rpmOrDiffuser_)
    {
        const vector avgCenter = gSum(p.Cf() * p.magSf())/totArea;
        const vector avgNormal = dir*gSum(p.Sf())/totArea;

        // Update angular velocity - convert [rpm] to [rad/s]
        tangentialVelocity =
            (rpm_*constant::mathematical::pi/30.0)
          * (p.Cf() - avgCenter) ^ avgNormal;

        normalVelocity = avgURemote;
    }
    else
    {
        normalVelocity = avgURemote;

        //revise direction
        nPF = diffuserDir_ / ((dir*p.nf()) & diffuserDir_);
    }

    if (minFlowRate_ >= 0.0)
    {
        normalVelocity = max(normalVelocity, minFlowRate_/totArea);
    }

    if (maxFlowRate_ > 0.0)
    {
        normalVelocity = min(normalVelocity, maxFlowRate_/totArea);
    }

    this->operator==(tangentialVelocity + nPF*normalVelocity);

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void pressureDirectedSwirlFvPatchField::write(Ostream& os) const
{
    // Note: do not write value
    fvPatchField<vector>::write(os);

    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);

    if (rpmOrDiffuser_)
    {
        writeEntry(os, "rpm", rpm_);
    }
    else
    {
        writeEntry(os, "diffuser", diffuserDir_);
    }

    if (minFlowRate_ >= 0.0)
    {
        writeEntry(os, "minFlowRate", minFlowRate_);
    }

    if (maxFlowRate_ > 0.0)
    {
        writeEntry(os, "maxFlowRate", maxFlowRate_);
    }

    writeEntry(os, "direction",
        ventilatorFlowDirectionNames_[direction_]
        );

    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
   fvPatchVectorField,
   pressureDirectedSwirlFvPatchField
);

} // End namespace Foam

// ************************************************************************* //
