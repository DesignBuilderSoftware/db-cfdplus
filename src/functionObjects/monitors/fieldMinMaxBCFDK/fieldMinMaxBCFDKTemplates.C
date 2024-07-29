/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's sampling and postProcessing libraries.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2013-2018 FSD blueCAPE Lda  http://www.bluecape.com.pt
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

#include "fieldMinMaxBCFDK.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


template<class Type>
void Foam::functionObjects::fieldMinMaxBCFDK::calcMinMaxFields
(
    const word& fieldName,
    const modeType& mode,
    const word& identifier
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    //true= min, false= max
    const bool minOrMax =
        this->operation_ == mdMinimum ? true : false;

    const bool checkInternalMesh =
        (this->dataMode_ == mdBoth) || (this->dataMode_ == mdCells);

    const bool checkPatches =
        (this->dataMode_ == mdBoth) || (this->dataMode_ == mdBoundaries);

    if (obr_.foundObject<fieldType>(fieldName))
    {
        const label proci = Pstream::myProcNo();

        const fieldType& origField = obr_.lookupObject<fieldType>(fieldName);
        const volScalarField field
        (
            convertField<volScalarField, fieldType>(origField, mode)
        );

        const volScalarField::Boundary& fieldBoundary =
            field.boundaryField();

        List<scalar> minVs(Pstream::nProcs());
        List<scalar> maxVs(Pstream::nProcs());

        if(checkInternalMesh)
        {
            if(minOrMax)
            {
                label minProcI = findMin(field);
                minVs[proci] = field[minProcI];
            }
            else
            {
                label maxProcI = findMax(field);
                maxVs[proci] = field[maxProcI];
            }
        }

        if(checkPatches)
        {
            forAll(fieldBoundary, patchI)
            {
                const scalarField& fp = fieldBoundary[patchI];
                if (fp.size())
                {
                    if(minOrMax)
                    {
                        label minPI = findMin(fp);
                        if (fp[minPI] < minVs[proci])
                        {
                            minVs[proci] = fp[minPI];
                        }
                    }
                    else
                    {
                        label maxPI = findMax(fp);
                        if (fp[maxPI] > maxVs[proci])
                        {
                            maxVs[proci] = fp[maxPI];
                        }
                    }
                }
            }
        }

        if(minOrMax)
        {
            Pstream::gatherList(minVs);
        }
        else
        {
            Pstream::gatherList(maxVs);
        }

        if (Pstream::master())
        {
            scalar result = 0.0;

            if(minOrMax)
            {
                label minI = findMin(minVs);
                result = minVs[minI];
            }
            else
            {
                label maxI = findMax(maxVs);
                result = maxVs[maxI];
            }

            Log
                << identifier
                << token::ASSIGN << token::SPACE
                << result << token::END_STATEMENT
                << token::SPACE;

            if (logToFile_)
            {
                file() << tab << result;
            }
        }
    }
}


// ************************************************************************* //
