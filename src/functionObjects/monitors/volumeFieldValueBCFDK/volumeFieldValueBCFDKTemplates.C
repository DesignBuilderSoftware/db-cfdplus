/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Derived from OpenFOAM's volFieldValue function object.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2013-2019 FSD blueCAPE Lda  http://www.bluecape.com.pt
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

#include "volumeFieldValueBCFDK.H"
#include "volFields.H"
#include "sampledSurface.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool
Foam::functionObjects::fieldValues::volumeFieldValueBCFDK::writeValuesBCFDK
(
    const word& fieldName,
    const modeType& mode,
    const word& identifier
)
{
    const bool ok = validField<Type>(fieldName);

    if (ok)
    {
        Field<Type> originalValues(setFieldValues<Type>(fieldName));
        scalarField values
        (
            convertField<scalarField,Field<Type>>(originalValues, mode)
        );

        scalarField V(filterField(fieldValue::mesh_.V()));
        scalarField weightField(values.size(), 1.0);

        forAll(weightFieldNames_, i)
        {
            weightField *= setFieldValues<scalar>(weightFieldNames_[i], true);
        }

        scalar result = processValues(values, V, weightField);

        if (Pstream::master())
        {
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

    return ok;
}

// ************************************************************************* //
