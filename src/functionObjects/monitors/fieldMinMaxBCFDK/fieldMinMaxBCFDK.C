/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's fieldMinMax function object.
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
#include "fieldTypes.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldMinMaxBCFDK, 0);
    addToRunTimeSelectionTable(functionObject, fieldMinMaxBCFDK, dictionary);

}
}

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::fieldMinMaxBCFDK::operationType,
    2
>::names[] =
{
    "maximum",
    "minimum"
};

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::fieldMinMaxBCFDK::dataModeType,
    3
>::names[] =
{
    "cells",
    "boundaries",
    "both"
};


const Foam::NamedEnum<Foam::functionObjects::fieldMinMaxBCFDK::operationType, 2>
Foam::functionObjects::fieldMinMaxBCFDK::operationTypeNames_;

const Foam::NamedEnum<Foam::functionObjects::fieldMinMaxBCFDK::dataModeType, 3>
Foam::functionObjects::fieldMinMaxBCFDK::dataModeTypeNames_;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::fieldMinMaxBCFDK::writeFileHeader
(
    const label i
)
{
    commonDictBCFDK::writeFileHeader
    (
        file(i),
        operationTypeNames_[operation_]
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldMinMaxBCFDK::fieldMinMaxBCFDK
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    commonDictBCFDK(dict),
    operation_(mdMaximum),
    dataMode_(mdBoth),
    fieldSet_()
{
    logFiles::resetName(typeName);

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldMinMaxBCFDK::~fieldMinMaxBCFDK()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldMinMaxBCFDK::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    commonDictBCFDK::read(dict);

    operation_ = operationTypeNames_
    [
        dict.lookupOrDefault<word>("operation", "maximum")
    ];

    dataMode_ = dataModeTypeNames_
    [
        dict.lookupOrDefault<word>("dataMode", "both")
    ];

    dict.lookup("fields") >> fieldSet_;

    return true;
}


bool Foam::functionObjects::fieldMinMaxBCFDK::execute()
{
    // Do nothing - only valid on write

    return true;
}


bool Foam::functionObjects::fieldMinMaxBCFDK::write()
{
    //Start the line tag
    Log << name() <<  token::COLON << token::SPACE;

    if (logToFile_ && Pstream::master())
    {
        logFiles::write();
        writeTime(file());
    }

    //Process all requests
    forAll(fieldSet_, fieldI)
    {
        calcMinMaxFields<scalar>(
            fieldSet_[fieldI],
            modes_[fieldI],
            identifiers_[fieldI]
            );
        calcMinMaxFields<vector>(
            fieldSet_[fieldI],
            modes_[fieldI],
            identifiers_[fieldI]
            );

        //These fields are not currently supported for bCFDK
        //calcMinMaxFields<sphericalTensor>(
        //    fieldSet_[fieldI],
        //    modes_[fieldI],
        //    identifiers_[fieldI]
        //    );
        //calcMinMaxFields<symmTensor>(
        //    fieldSet_[fieldI],
        //    modes_[fieldI],
        //    identifiers_[fieldI]
        //    );
        //calcMinMaxFields<tensor>(
        //    fieldSet_[fieldI],
        //    modes_[fieldI],
        //    identifiers_[fieldI]
        //    );
    }

    if (logToFile_ && Pstream::master())
    {
        file() << endl;
    }

    //End line twice
    Log << nl << endl;

    return true;
}



// ************************************************************************* //
