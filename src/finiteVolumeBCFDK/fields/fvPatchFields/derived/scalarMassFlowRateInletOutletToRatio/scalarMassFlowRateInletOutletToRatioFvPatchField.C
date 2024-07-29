/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's fixedProfile.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2019-2020 FSD blueCAPE Lda  http://www.bluecape.com.pt
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

#include "scalarMassFlowRateInletOutletToRatioFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    // declare specialization within 'Foam' namespace
    template<>
    const char* NamedEnum
    <
        Foam::scalarMassFlowRateInletOutletToRatioFvPatchField::profileTypes,
        2
    >::names[] =
    {
        "massFlowRate",
        "fixedMassFraction"
    };

    template<>
    const char* NamedEnum
    <
        Foam::scalarMassFlowRateInletOutletToRatioFvPatchField
            ::massFractionTypes,
        2
    >::names[] =
    {
        "massFraction",
        "massRatio"
    };

    // * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * //

    const NamedEnum
    <
        scalarMassFlowRateInletOutletToRatioFvPatchField::profileTypes,
        2
    > scalarMassFlowRateInletOutletToRatioFvPatchField::profileTypeNames_;

    const NamedEnum
    <
        scalarMassFlowRateInletOutletToRatioFvPatchField::massFractionTypes,
        2
    > scalarMassFlowRateInletOutletToRatioFvPatchField::massFractionTypeNames_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::scalarMassFlowRateInletOutletToRatioFvPatchField
    ::scalarMassFlowRateInletOutletToRatioFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    profile_(),
    origin_(0),
    direction_(Zero),
    profileType_(fixedMassFraction),
    massFractionType_(massFraction),
    phiName_("phi"),
    rhoName_("rhoInf"),
    inletOnly_(false)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


Foam::scalarMassFlowRateInletOutletToRatioFvPatchField
    ::scalarMassFlowRateInletOutletToRatioFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
    profile_(Function1<scalar>::New("profile", dict)),
    origin_(dict.lookup<scalar>("origin")),
    direction_(dict.lookup("direction")),
    profileType_(profileTypeNames_.read(dict.lookup("profileType"))),
    massFractionType_
    (
        massFractionTypeNames_.read(dict.lookup("massFractionType"))
    ),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rhoInf")),
    inletOnly_(dict.lookupOrDefault<Switch>("inletOnly", false))
{
    if (mag(direction_) < SMALL)
    {
        FatalErrorInFunction
            << "The direction must be non-zero"
            << abort(FatalError);
    }

    // Normalise the direction
    direction_ /= mag(direction_);

    //Initialize the gradient
    this->refGrad() = 0.0;

    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );

        this->refValue() = *this;
    }
    else
    {
        // Evaluate profile
        evaluate(Pstream::commsTypes::blocking);
    }
}


Foam::scalarMassFlowRateInletOutletToRatioFvPatchField
    ::scalarMassFlowRateInletOutletToRatioFvPatchField
(
    const scalarMassFlowRateInletOutletToRatioFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(ptf, p, iF, mapper), //Have to map
    profile_(ptf.profile_, false),
    origin_(ptf.origin_),
    direction_(ptf.direction_),
    profileType_(ptf.profileType_),
    massFractionType_(ptf.massFractionType_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    inletOnly_(ptf.inletOnly_)
{
    // Cannot evaluate profile, since it depends on phi
    //this->evaluate();
}


Foam::scalarMassFlowRateInletOutletToRatioFvPatchField
    ::scalarMassFlowRateInletOutletToRatioFvPatchField
(
    const scalarMassFlowRateInletOutletToRatioFvPatchField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(ptf, iF),
    profile_(ptf.profile_, false),
    origin_(ptf.origin_),
    direction_(ptf.direction_),
    profileType_(ptf.profileType_),
    massFractionType_(ptf.massFractionType_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    inletOnly_(ptf.inletOnly_)
{
    // Evaluate the profile if defined... er, can't do it, because
    // reconstruction doesn't work with this, due to depending on phi
    //if (ptf.profile_.valid())
    //{
    //    this->evaluate();
    //}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::scalarMassFlowRateInletOutletToRatioFvPatchField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const scalarField dirCmpt
    (
        (direction_ & this->patch().Cf())
      - origin_
    );

    //Get velocity field time area on patch...
    //might as well get the field, so that it can also be used for knowing
    //which flow direction to follow
    const fvsPatchScalarField& phip = patch().
        lookupPatchField<surfaceScalarField, scalar>(this->phiName_);

    switch (profileType_)
    {
        case massFlowRate:
        {
            //Calculate the fluid's mass flow rate

            //Need switch negative sign, given that the flow rate is negative
            //by default when going into the domain
            fvsPatchScalarField fluidMassFlowRate(phip);

            fluidMassFlowRate *= scalar(-1);

            //Get density
            if (rhoName_=="rhoInf")
            {
                const IOdictionary& transportProperties =
                    db().lookupObject<IOdictionary>("transportProperties");
                const dimensionedScalar rho
                (
                    "rho",
                    dimDensity,
                    transportProperties
                );

                fluidMassFlowRate *= rho.value();
            }
            else if (rhoName_!="none")
            {
                const scalarField& prho =
                    patch().lookupPatchField<volScalarField, scalar>(rhoName_);

                fluidMassFlowRate *= prho;
            }

            //calculate the area weighted profile of the contaminant mass flow
            //rate
            scalarField scalarMassFlowRate
            (
                profile_->value(dirCmpt)
              * this->patch().magSf() / gSum(this->patch().magSf())
            );

            switch (massFractionType_)
            {
                case massFraction:
                {
                    this->refValue() =
                    (
                        scalarMassFlowRate
                      / (scalarMassFlowRate + fluidMassFlowRate)
                    );
                    break;
                }
                case massRatio:
                {
                    this->refValue() =
                    (
                        scalarMassFlowRate
                      / fluidMassFlowRate
                    );
                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Unknown mass fraction type. Valid types are: "
                        << massFractionTypeNames_ << nl
                        << exit(FatalError);
                }
            }

            break;
        }
        case fixedMassFraction:
        {
            //Use the profile as-is
            //Does not rely on mass fraction type
            this->refValue() = profile_->value(dirCmpt);
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown profile type. Valid types are: "
                << profileTypeNames_ << nl
                << exit(FatalError);
        }
    }

    if(inletOnly_)
    {
        this->valueFraction() = 1.0;
    }
    else
    {
        this->valueFraction() = 1.0 - pos0(phip);
    }

    inletOutletFvPatchScalarField::updateCoeffs();
}


void Foam::scalarMassFlowRateInletOutletToRatioFvPatchField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, profile_());
    writeEntry(os, "direction", direction_);
    writeEntry(os, "origin", origin_);
    writeEntry(os, "profileType", profileTypeNames_[profileType_]);
    writeEntry(os, "massFractionType",
       massFractionTypeNames_[massFractionType_]
       );

    if (phiName_ != "phi")
    {
        writeEntry(os, "phi", phiName_);
    }

    if (rhoName_ != "rhoInf")
    {
        writeEntry(os, "rho", rhoName_);
    }

    writeEntry(os, "inletOnly", inletOnly_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        scalarMassFlowRateInletOutletToRatioFvPatchField
    );
}

// ************************************************************************* //
