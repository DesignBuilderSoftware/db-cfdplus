/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "wallHeatFluxTemperatureBCFDKFvPatchScalarField.H"
#include "volFields.H"
#include "physicoChemicalConstants.H"
#include "addToRunTimeSelectionTable.H"

using Foam::constant::physicoChemical::sigma;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char*
    NamedEnum
    <
        wallHeatFluxTemperatureBCFDKFvPatchScalarField::operationMode,
        2
    >::names[] =
    {
        "power",
        "flux"
    };

    template<>
    const char*
    NamedEnum
    <
        wallHeatFluxTemperatureBCFDKFvPatchScalarField::operationHtcMode,
        3
    >::names[] =
    {
        "zeroGradient",
        "defined",
        "fluidKappa"
    };
}

const Foam::NamedEnum
<
    Foam::wallHeatFluxTemperatureBCFDKFvPatchScalarField::operationMode,
    2
> Foam::wallHeatFluxTemperatureBCFDKFvPatchScalarField::operationModeNames;

const Foam::NamedEnum
<
    Foam::wallHeatFluxTemperatureBCFDKFvPatchScalarField::operationHtcMode,
    3
> Foam::wallHeatFluxTemperatureBCFDKFvPatchScalarField::htcModeNames;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallHeatFluxTemperatureBCFDKFvPatchScalarField::
wallHeatFluxTemperatureBCFDKFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch()),
    mode_(fixedHeatFlux),
    Q_(0),
    q_(),
    h_(0),
    Tmin_(0),
    Tmax_(great)
{
}


Foam::wallHeatFluxTemperatureBCFDKFvPatchScalarField::
wallHeatFluxTemperatureBCFDKFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    mode_(operationModeNames.read(dict.lookup("mode"))),
    htcMode_(htcModeNames.read(dict.lookup("htcMode"))),
    Q_(0),
    q_(),
    h_(0),
    Tmin_(0),
    Tmax_(great)
{
    switch (mode_)
    {
        case fixedPower:
        {
            dict.lookup("Q") >> Q_;

            break;
        }
        case fixedHeatFlux:
        {
            q_ = scalarField("q", dict, p.size());

            break;
        }
    }

    switch (htcMode_)
    {
        case zeroGradient:
        {
            break;
        }
        case defined:
        {
            dict.lookup("h") >> h_;

            break;
        }
        case fluidKappa:
        {
            break;
        }
    }

    Tmin_ = dict.lookupOrDefault("Tmin", 0.0);
    Tmax_ = dict.lookupOrDefault("Tmax", great);

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
}


Foam::wallHeatFluxTemperatureBCFDKFvPatchScalarField::
wallHeatFluxTemperatureBCFDKFvPatchScalarField
(
    const wallHeatFluxTemperatureBCFDKFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    temperatureCoupledBase(patch(), ptf),
    mode_(ptf.mode_),
    htcMode_(ptf.htcMode_),
    Q_(ptf.Q_),
    q_(),
    h_(ptf.h_),
    Tmin_(ptf.Tmin_),
    Tmax_(ptf.Tmax_)
{
    switch (mode_)
    {
        case fixedPower:
        {
            break;
        }
        case fixedHeatFlux:
        {
            mapper(q_, ptf.q_);
            break;
        }
    }
}


Foam::wallHeatFluxTemperatureBCFDKFvPatchScalarField::
wallHeatFluxTemperatureBCFDKFvPatchScalarField
(
    const wallHeatFluxTemperatureBCFDKFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    temperatureCoupledBase(tppsf),
    mode_(tppsf.mode_),
    htcMode_(tppsf.htcMode_),
    Q_(tppsf.Q_),
    q_(tppsf.q_),
    h_(tppsf.h_),
    Tmin_(tppsf.Tmin_),
    Tmax_(tppsf.Tmax_)
{}


Foam::wallHeatFluxTemperatureBCFDKFvPatchScalarField::
wallHeatFluxTemperatureBCFDKFvPatchScalarField
(
    const wallHeatFluxTemperatureBCFDKFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    temperatureCoupledBase(patch(), tppsf),
    mode_(tppsf.mode_),
    htcMode_(tppsf.htcMode_),
    Q_(tppsf.Q_),
    q_(tppsf.q_),
    h_(tppsf.h_),
    Tmin_(tppsf.Tmin_),
    Tmax_(tppsf.Tmax_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wallHeatFluxTemperatureBCFDKFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);

    switch (mode_)
    {
        case fixedPower:
        {
            break;
        }
        case fixedHeatFlux:
        {
            m(q_, q_);

            break;
        }
    }
}


void Foam::wallHeatFluxTemperatureBCFDKFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const wallHeatFluxTemperatureBCFDKFvPatchScalarField& tiptf =
        refCast<const wallHeatFluxTemperatureBCFDKFvPatchScalarField>(ptf);

    switch (mode_)
    {
        case fixedPower:
        {
            break;
        }
        case fixedHeatFlux:
        {
            q_.rmap(tiptf.q_, addr);

            break;
        }
    }
}


void Foam::wallHeatFluxTemperatureBCFDKFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalarField& Tp = *this;

    switch (mode_)
    {
        case fixedPower:
        {
            q_ = Q_/gSum(patch().magSf());

            break;
        }
        case fixedHeatFlux:
        {
            break;
        }
    }

    const scalarField Tc(this->patchInternalField());

    switch (htcMode_)
    {
        case zeroGradient:
        {
            Tp = Tc;
            break;
        }
        case defined:
        {
            Tp = (q_ + h_*Tc) / h_;

            break;
        }
        case fluidKappa:
        {
            const scalarField kappaDeltaCoeffs
            (
                this->kappa(*this)*patch().deltaCoeffs()
            );

            Info << "kappaDeltaCoeffs: " << kappaDeltaCoeffs << endl;

            Tp = (q_ + kappaDeltaCoeffs*Tc) / kappaDeltaCoeffs;

            break;
        }
    }

    Tp = max(min(Tp, Tmax_), Tmin_);

    fixedValueFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        const scalar Q = gSum(kappa(*this)*patch().magSf()*snGrad());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " :"
            << " estimated heat transfer rate:" << Q
            << ", wall temperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }
}


void Foam::wallHeatFluxTemperatureBCFDKFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);

    writeEntry(os, "mode", operationModeNames[mode_]);
    temperatureCoupledBase::write(os);

    switch (mode_)
    {
        case fixedPower:
        {
            writeEntry(os, "Q", Q_);

            break;
        }
        case fixedHeatFlux:
        {
            writeEntry(os, "q", q_);

            break;
        }
    }

    writeEntry(os, "htcMode", htcModeNames[htcMode_]);

    switch (htcMode_)
    {
        case zeroGradient:
        {
            break;
        }
        case defined:
        {
            writeEntry(os, "h", h_);

            break;
        }
        case fluidKappa:
        {
            break;
        }
    }

    writeEntryIfDifferent(os, "Tmin", scalar(0), Tmin_);
    writeEntryIfDifferent(os, "Tmax", great, Tmax_);

    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        wallHeatFluxTemperatureBCFDKFvPatchScalarField
    );
}

// ************************************************************************* //
