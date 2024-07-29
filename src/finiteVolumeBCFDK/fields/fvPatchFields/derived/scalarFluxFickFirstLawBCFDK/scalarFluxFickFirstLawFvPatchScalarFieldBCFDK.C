/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Somewhat based on OpenFOAM's old turbulentHeatFluxTemperature, now named
        externalWallHeatFluxTemperature.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2014-2021 FSD blueCAPE Lda  http://www.bluecape.com.pt
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

#include "scalarFluxFickFirstLawFvPatchScalarFieldBCFDK.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    // declare specialization within 'Foam' namespace
    template<>
    const char* NamedEnum
    <
        Foam::scalarFluxFickFirstLawFvPatchScalarFieldBCFDK::massFractionTypes,
        2
    >::names[] =
    {
        "massFraction",
        "massRatio"
    };

    // * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * //

    const NamedEnum
    <
        scalarFluxFickFirstLawFvPatchScalarFieldBCFDK::massFractionTypes,
        2
    > scalarFluxFickFirstLawFvPatchScalarFieldBCFDK::massFractionTypeNames_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::scalarFluxFickFirstLawFvPatchScalarFieldBCFDK::
scalarFluxFickFirstLawFvPatchScalarFieldBCFDK
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    fluxValue_(p.size(), 0.0),
    molDiffEffName_("undefinedMolDiffEff"),
    rhoName_("rhoInf"),
    massFractionType_(massFraction)
{}


Foam::scalarFluxFickFirstLawFvPatchScalarFieldBCFDK::
scalarFluxFickFirstLawFvPatchScalarFieldBCFDK
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    fluxValue_("flux", dict, p.size()),
    molDiffEffName_(dict.lookup("molDiffEff")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rhoInf")),
    massFractionType_
    (
        massFractionTypeNames_.read(dict.lookup("massFractionType"))
    )
{
    fvPatchScalarField::operator=(patchInternalField());

    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


Foam::scalarFluxFickFirstLawFvPatchScalarFieldBCFDK::
scalarFluxFickFirstLawFvPatchScalarFieldBCFDK
(
    const scalarFluxFickFirstLawFvPatchScalarFieldBCFDK& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    fluxValue_(mapper(ptf.fluxValue_)),
    molDiffEffName_(ptf.molDiffEffName_),
    rhoName_(ptf.rhoName_),
    massFractionType_(ptf.massFractionType_)
{}


Foam::scalarFluxFickFirstLawFvPatchScalarFieldBCFDK::
scalarFluxFickFirstLawFvPatchScalarFieldBCFDK
(
    const scalarFluxFickFirstLawFvPatchScalarFieldBCFDK& thftpsf
)
:
    mixedFvPatchScalarField(thftpsf),
    fluxValue_(thftpsf.fluxValue_),
    molDiffEffName_(thftpsf.molDiffEffName_),
    rhoName_(thftpsf.rhoName_),
    massFractionType_(thftpsf.massFractionType_)
{}


Foam::scalarFluxFickFirstLawFvPatchScalarFieldBCFDK::
scalarFluxFickFirstLawFvPatchScalarFieldBCFDK
(
    const scalarFluxFickFirstLawFvPatchScalarFieldBCFDK& thftpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(thftpsf, iF),
    fluxValue_(thftpsf.fluxValue_),
    molDiffEffName_(thftpsf.molDiffEffName_),
    rhoName_(thftpsf.rhoName_),
    massFractionType_(thftpsf.massFractionType_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::scalarFluxFickFirstLawFvPatchScalarFieldBCFDK::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
    m(fluxValue_, fluxValue_);
}


void Foam::scalarFluxFickFirstLawFvPatchScalarFieldBCFDK::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const scalarFluxFickFirstLawFvPatchScalarFieldBCFDK& thftptf =
        refCast<const scalarFluxFickFirstLawFvPatchScalarFieldBCFDK>
        (
            ptf
        );

    fluxValue_.rmap(thftptf.fluxValue_, addr);
}


void Foam::scalarFluxFickFirstLawFvPatchScalarFieldBCFDK::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField& molDiffEffp =
        patch().lookupPatchField<volScalarField, scalar>
        (
            molDiffEffName_
        );

    scalarField rho(patch().size());

    //Get density
    if (rhoName_=="rhoInf")
    {
        const IOdictionary& transportProperties =
            db().lookupObject<IOdictionary>("transportProperties");
        const dimensionedScalar rhoInf
        (
            "rho",
            dimDensity,
            transportProperties
        );

        rho = rhoInf.value();
    }
    else
    {
        rho = patch().lookupPatchField<volScalarField, scalar>(rhoName_);
    }

    tmp<scalarField> ttransportedField = this->patchInternalField();
    scalarField& transportedField = ttransportedField.ref();

    switch (massFractionType_)
    {
        case massFraction:
        {
            //Need to calculate the total mass density, while taking into
            //account the density of the scalar, in function of the mass
            //fraction on the adjacent cells

            rho += rho*transportedField/(scalar(1)-transportedField);
            break;
        }
        case massRatio:
        {
            //No change required, the base fluid density is enough
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

    refGrad() = fluxValue_ / (molDiffEffp*rho);

    this->valueFraction() =
        1.0 - pos0(transportedField + refGrad()/this->patch().deltaCoeffs());

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::scalarFluxFickFirstLawFvPatchScalarFieldBCFDK::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);

    writeEntry(os, "flux", fluxValue_);
    writeEntry(os, "molDiffEff", molDiffEffName_);

    if (rhoName_ != "rhoInf")
    {
        writeEntry(os, "rho", rhoName_);
    }

    writeEntry(os, "massFractionType",
        massFractionTypeNames_[massFractionType_]
        );

    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        scalarFluxFickFirstLawFvPatchScalarFieldBCFDK
    );
}

// ************************************************************************* //

