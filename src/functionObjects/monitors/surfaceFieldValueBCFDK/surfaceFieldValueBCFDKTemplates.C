/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Derived from OpenFOAM's surfaceFieldValue function object.
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

#include "surfaceFieldValueBCFDK.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "sampledSurface.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool
Foam::functionObjects::fieldValues::surfaceFieldValueBCFDK::writeValuesBCFDK
(
    const word& fieldName,
    const scalarField& weightField,
    const bool orient,
    const modeType& mode,
    const word& identifier
)
{
    const bool ok = validField<Type>(fieldName);

    if (ok)
    {
        Field<Type> originalValues
        (
            getFieldValues<Type>(fieldName, true, orient)
        );
        scalarField values
        (
            convertField<scalarField,Field<Type>>(originalValues, mode)
        );

        vectorField Sf;
        if (surfacePtr_.valid())
        {
            // Get oriented Sf
            Sf = surfacePtr_().Sf();
        }
        else
        {
            // Get oriented Sf
            Sf = filterField(mesh_.Sf(), true);
        }

        // Combine onto master
        combineFields(values);
        combineFields(Sf);

        if (operation_ != operationType::none)
        {
            // Apply scale factor
            values *= scaleFactor_;

            if (Pstream::master())
            {
                scalar result = processValues(values, Sf, weightField);

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

    return ok;
}

// ************************************************************************* //
