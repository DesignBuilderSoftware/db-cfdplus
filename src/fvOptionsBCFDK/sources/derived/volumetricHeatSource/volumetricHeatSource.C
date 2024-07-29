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

#include "volumetricHeatSource.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    template<>
    const char* NamedEnum<fv::volumetricHeatSource::heatSourceType, 2>::
    names[] =
    {
        "power",
        "rate"
    };

    namespace fv
    {
        defineTypeNameAndDebug(volumetricHeatSource, 0);
        addToRunTimeSelectionTable
        (
            option,
            volumetricHeatSource,
            dictionary
        );

        // * * * * * * * * * * * * Static Data Members * * * * * * * * * * * //

        const NamedEnum
        <
            volumetricHeatSource::heatSourceType,
            2
        > volumetricHeatSource::heatSourceTypeNames_;
    }
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::volumetricHeatSource::volumetricHeatSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    heatSource_(heatSourceTypeNames_.read(coeffs_.lookup("heatSource"))),
    q_(NULL)
{
    switch (heatSource_)
    {
        case hsPower:
        case hsFlux:
        {
            q_.reset
            (
                Function1<scalar>::New("q", coeffs_).ptr()
            );
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown heat source type. Valid types are: "
                << heatSourceTypeNames_ << nl
                << exit(FatalError);
        }
    }

    const basicThermo& thermo =
        mesh_.lookupObject<basicThermo>(basicThermo::dictName);

    fieldNames_.setSize(1, thermo.he().name());
    applied_.setSize(1, false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::volumetricHeatSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label
)
{
    const scalarField& V = mesh_.V();
    const scalar t = mesh_.time().value();
    scalar q_final(q_->value(t));
    scalarField hp_source(cells_.size());

    if (this->V() > VSMALL && mag(q_final) > VSMALL)
    {
        // Note: hp_source is received by OpenFOAM and then it divides by
        // the volume of each cell. That is why:
        //  - With heat flux, we must negate the cell volume division
        //    that OpenFOAM does.
        //  - When using power, we must first convert to flux and then can
        //    calculate the same way as with flux.

        switch (heatSource_)
        {
            case hsPower:
            {
                forAll(cells_, i)
                {
                    label cellI = cells_[i];
                    hp_source[i] =
                        (
                            q_final * V[cellI]
                        )
                      /
                        this->V();
                }
                break;
            }
            case hsFlux:
            {
                forAll(cells_, i)
                {
                    label cellI = cells_[i];
                    hp_source[i] =
                        (
                            q_final * V[cellI]
                        );
                }
                break;
            }
            default:
            {}
        }

        scalarField& heSource = eqn.source();

        if (debug)
        {
            Info << "Source term: " << heSource << endl;
        }

        forAll(cells_, i)
        {
            heSource[cells_[i]] -= hp_source[i];
        }
    }
}


bool Foam::fv::volumetricHeatSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        if (coeffs_.found(q_->name()))
        {
            q_.reset
            (
                Function1<scalar>::New(q_->name(), dict).ptr()
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
