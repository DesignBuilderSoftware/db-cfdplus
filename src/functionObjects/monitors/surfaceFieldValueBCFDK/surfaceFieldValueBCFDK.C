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
#include "fvMesh.H"
#include "cyclicPolyPatch.H"
#include "emptyPolyPatch.H"
#include "coupledPolyPatch.H"
#include "sampledSurface.H"
#include "mergePoints.H"
#include "indirectPrimitivePatch.H"
#include "addToRunTimeSelectionTable.H"

#if defined( WIN32 ) || defined( WIN64 )
#include "PatchTools.T.H"
#else
#include "PatchTools.H"
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace fieldValues
{
    defineTypeNameAndDebug(surfaceFieldValueBCFDK, 0);
    addToRunTimeSelectionTable(fieldValue, surfaceFieldValueBCFDK, dictionary);
    addToRunTimeSelectionTable
    (
        functionObject,
        surfaceFieldValueBCFDK,
        dictionary
    );
}
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::surfaceFieldValueBCFDK::
surfaceFieldValueBCFDK
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    Foam::functionObjects::fieldValues::surfaceFieldValue
    (
        name,
        runTime,
        dict
    ),
    commonDictBCFDK(dict)
{
    read(dict);
}

Foam::functionObjects::fieldValues::surfaceFieldValueBCFDK::
surfaceFieldValueBCFDK
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    Foam::functionObjects::fieldValues::surfaceFieldValue
    (
        name,
        obr,
        dict
    ),
    Foam::commonDictBCFDK::commonDictBCFDK(dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::surfaceFieldValueBCFDK::
~surfaceFieldValueBCFDK()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldValues::surfaceFieldValueBCFDK::read
(
    const dictionary& dict
)
{
    surfaceFieldValue::read(dict);

    commonDictBCFDK::read(dict);

    return true;
}


bool Foam::functionObjects::fieldValues::surfaceFieldValueBCFDK::write()
{
    //Note: We removed all writing to disk, since we don't use it.

    //Start the line tag
    Log << name() <<  token::COLON << token::SPACE;

    if (logToFile_ && Pstream::master())
    {
        logFiles::write();
        writeTime(file());
    }

    // Construct weight field. Note: zero size means weight = 1
    scalarField weightField;
    if (weightFieldName_ != "none")
    {
        weightField =
            getFieldValues<scalar>
            (
                weightFieldName_,
                true,
                orientWeightField_
            );
    }

    // Combine onto master
    combineFields(weightField);

    // Process the fields
    forAll(fields_, i)
    {
        const word& fieldName = fields_[i];
        const modeType& modeName = modes_[i];
        const word& idName = identifiers_[i];
        bool ok = false;

        bool orient = i >= orientedFieldsStart_;
        ok = ok || writeValuesBCFDK<scalar>
            (
                fieldName,
                weightField,
                orient,
                modeName,
                idName
            );
        ok = ok || writeValuesBCFDK<vector>
            (
                fieldName,
                weightField,
                orient,
                modeName,
                idName
            );
        ok = ok || writeValuesBCFDK<sphericalTensor>
            (
                fieldName,
                weightField,
                orient,
                modeName,
                idName
            );
        ok = ok || writeValuesBCFDK<symmTensor>
            (
                fieldName,
                weightField,
                orient,
                modeName,
                idName
            );
        ok = ok || writeValuesBCFDK<tensor>
            (
                fieldName,
                weightField,
                orient,
                modeName,
                idName
            );

        if (!ok)
        {
            WarningInFunction
                << "Requested field " << fieldName
                << " not found in database and not processed"
                << endl;
        }
    }

    if (logToFile_ && Pstream::master())
    {
        file() << endl;
    }

    Log << nl << endl;

    return true;
}


// ************************************************************************* //
