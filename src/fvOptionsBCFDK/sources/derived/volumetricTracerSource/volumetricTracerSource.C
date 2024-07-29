/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's effectivenessHeatExchangerSource.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2014-2019 FSD blueCAPE Lda  http://www.bluecape.com.pt
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

#include "volumetricTracerSource.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(volumetricTracerSource, 0);
        addToRunTimeSelectionTable
        (
            option,
            volumetricTracerSource,
            dictionary
        );

        template<>
        const char* NamedEnum
        <
            volumetricTracerSource::sourceTypes,
            3
        >::names[] =
        {
            "constantMassFraction",
            "massFlowRate",
            "massFractionFlowRate"
        };

        template<>
        const char* NamedEnum
        <
            volumetricTracerSource::massFractionTypes,
            2
        >::names[] =
        {
            "massFraction",
            "massRatio"
        };

    }
}


const Foam::NamedEnum
<
    Foam::fv::volumetricTracerSource::sourceTypes,
    3
> Foam::fv::volumetricTracerSource::sourceTypeNames_;

const Foam::NamedEnum
<
    Foam::fv::volumetricTracerSource::massFractionTypes,
    2
> Foam::fv::volumetricTracerSource::massFractionTypeNames_;


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::volumetricTracerSource::volumetricTracerSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    tSource_(NULL),
    tracerName_(coeffs_.lookup("tracerName")),
    sourceType_(sourceTypeNames_.read(coeffs_.lookup("sourceType"))),
    rhoName_(coeffs_.lookup("rho")),
    massFractionType_
    (
        massFractionTypeNames_.read(coeffs_.lookup("massFractionType"))
    )
{
    tSource_.reset
    (
        Function1<scalar>::New("tSource", coeffs_).ptr()
    );

    fieldNames_.setSize(1, tracerName_);
    applied_.setSize(1, false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::volumetricTracerSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label
)
{
    //Don't calculate here those that are for constrain
    if(sourceType_ == constantMassFraction) return;

    const scalar t = mesh_.time().value();
    scalar tSource_final(tSource_->value(t));

    // The tracer source needs to be in tracer's units per second

    if (this->V() > VSMALL && mag(tSource_final) > VSMALL)
    {
        scalarField& tracerSource = eqn.source();
        const scalarField& V = mesh_.V();

        scalarField rhoFinal(cells_.size(), scalar(1.0));
        if (rhoName_=="rhoInf")
        {
            const dictionary& transportProperties =
                mesh_.lookupObject<dictionary>("transportProperties");
            const dimensionedScalar rhoInf
            (
                "rho",
                dimDensity,
                transportProperties
            );

            rhoFinal = rhoInf.value();
        }
        else if (rhoName_!="none")
        {
            const volScalarField& rho =
                mesh_.lookupObject<volScalarField>(rhoName_);

            forAll(cells_, i)
            {
                label cellI = cells_[i];
                rhoFinal[i] = rho[cellI];
            }
        }

        if(sourceType_ == massFractionFlowRate)
        {
            // mass fraction per volume per second, needs us to first multiply
            // by the total volume of the cell zone
            tSource_final *= this->V();

            // Source term needs to be multiplied by the volume of the cell
            // because it is being added to the equation matrix, which also
            // takes into account the volume of the cell.
            forAll(cells_, i)
            {
                label cellI = cells_[i];
                tracerSource[cellI] -= tSource_final*V[cellI]*rhoFinal[i];
            }
        }
        else if(sourceType_ == massFlowRate)
        {

            // Source term is divided by the density, so it depends on the
            // origin for the base fluid density and the mass fraction type
            // In other words: if the equation is divided by rho, then we must
            // also divide the source term by rho...

            switch (massFractionType_)
            {
                case massFraction:
                {
                    //Need to calculate the total mass density, while taking
                    //into account the density of the scalar, in function of
                    //the mass fraction on the adjacent cells
                    const volScalarField& transportedField = eqn.psi();

                    forAll(cells_, i)
                    {
                        label cellI = cells_[i];

                        rhoFinal[i] +=
                            rhoFinal[i]
                          * (
                                transportedField[cellI]
                              / (scalar(1) - transportedField[cellI])
                            );
                    }

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


            // Source term needs to be multiplied by the volume of the cell
            // because it is being added to the equation matrix, which also
            // takes into account the volume of the cell.
            forAll(cells_, i)
            {
                label cellI = cells_[i];
                tracerSource[cellI] -= tSource_final*V[cellI] / rhoFinal[i];
            }

        }
    }
}


void Foam::fv::volumetricTracerSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    // Simply call the non-rho implementation, which will go get the rho field
    // we want
    addSup(eqn, fieldi);
}


void Foam::fv::volumetricTracerSource::constrain
(
    fvMatrix<scalar>& eqn,
    const label
)
{
    //Don't calculate here those that are not for constrain
    if(sourceType_ != constantMassFraction) return;

    DebugInfo
        << "volumetricTracerSource::constrain for source "
        << tracerName_
        << endl;

    const scalar t = mesh_.time().value();

    // The tracer source needs to be in tracer's units per second

    if (this->V() > VSMALL)
    {
        scalarField tSource_final(cells_.size(), tSource_->value(t));

        //This is because the value comes in as constant mass fraction per unit
        //volume
        tSource_final *= this->V();

        eqn.setValues(cells_, tSource_final);
    }
}


bool Foam::fv::volumetricTracerSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        if (coeffs_.found(tSource_->name()))
        {
            tSource_.reset
            (
                Function1<scalar>::New(tSource_->name(), dict).ptr()
            );
        }

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
