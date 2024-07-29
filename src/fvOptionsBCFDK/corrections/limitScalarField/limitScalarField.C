/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's limitVelocity.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2021 FSD blueCAPE Lda  http://www.bluecape.com.pt
-------------------------------------------------------------------------------
License
    This file is a derivative of OpenFOAM.

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

#include "limitScalarField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(limitScalarField, 0);
    addToRunTimeSelectionTable
    (
        option,
        limitScalarField,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::limitScalarField::limitScalarField
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    targetFieldName_(coeffs_.lookup<word>("targetField")),
    refFieldName_(coeffs_.lookup<word>("refField")),
    scaleMax_(coeffs_.lookup<scalar>("scaleMax"))
{
    fieldNames_.setSize(1, targetFieldName_);
    applied_.setSize(1, false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::limitScalarField::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        coeffs_.lookup("scaleMax") >> scaleMax_;

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::fv::limitScalarField::correct(volScalarField& field)
{
    const volScalarField& refField =
        mesh_.lookupObject<volScalarField>(refFieldName_);

    const scalar maxValue = max(refField.primitiveField())*scaleMax_;

    scalarField& targetFr = field.primitiveFieldRef();

    forAll(cells_, i)
    {
        const label celli = cells_[i];

        if (targetFr[celli] > maxValue)
        {
            targetFr[celli] = maxValue;
        }
    }

    // handle boundaries in the case of 'all'
    if (selectionMode_ == smAll)
    {
        volScalarField::Boundary& targetBFr = field.boundaryFieldRef();

        forAll(targetBFr, patchi)
        {
            fvPatchScalarField& patchSF = targetBFr[patchi];

            if (!patchSF.fixesValue())
            {
                forAll(patchSF, facei)
                {
                    if (patchSF[facei] > maxValue)
                    {
                        patchSF[facei] = maxValue;
                    }
                }
            }
        }
    }
}


// ************************************************************************* //
