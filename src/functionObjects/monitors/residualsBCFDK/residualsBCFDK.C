/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Derived from OpenFOAM's residuals function object.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2023-2023 FSD blueCAPE Lda  http://www.bluecape.com.pt
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

#include "residualsBCFDK.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(residualsBCFDK, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        residualsBCFDK,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::residualsBCFDK::residualsBCFDK
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    fieldSet_(),
    fieldSetSubStep_(),
    writeToFile_(true)
{
    read(dict);

    if (writeToFile_)
    {
        // Compile list of file names:
        //  - first absolute residuals
        //  - second original ones

        wordList residualFileNames;
        word origResidualName(typeName);
        origResidualName.replace("BCFDK","");

        residualFileNames.append(typeName);
        residualFileNames.append(origResidualName);

        resetNames(residualFileNames);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::residualsBCFDK::~residualsBCFDK()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::residualsBCFDK::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    debug = dict.lookupOrDefault<label>("debug", debug);

    writeToFile_ = dict.lookupOrDefault<bool>("writeToFile", true);

    dict.lookup("fields") >> fieldSet_;

    fieldSetSubStep_ = dict.lookupOrDefault("substeps", labelList());

    fieldTranslationToMatrixSolver_  =
        dict.lookupOrDefault("fieldToSolver", HashTable<word>());

    if (fieldSetSubStep_.size() == 0)
    {
        //reset to be all zeros
        fieldSetSubStep_.setSize(fieldSet_.size(), 0);
    }
    else if (fieldSetSubStep_.size() != fieldSet_.size())
    {
        FatalErrorInFunction
            << "The list of substeps does not have the same size as the list"
               " of fields."
            << exit(FatalError);
    }

    return true;
}


void Foam::functionObjects::residualsBCFDK::writeFileHeader(const label i)
{
    if (Pstream::master() && writeToFile_)
    {
        writeHeader(file(i), "residualsBCFDK");
        writeCommented(file(i), "Time");

        forAll(fieldSet_, fieldi)
        {
            const word& fieldName = fieldSet_[fieldi];
            const label substep = fieldSetSubStep_[fieldi];

            writeFileHeader<scalar>(i, fieldName, substep);
            writeFileHeader<vector>(i, fieldName, substep);
            writeFileHeader<sphericalTensor>(i, fieldName, substep);
            writeFileHeader<symmTensor>(i, fieldName, substep);
            writeFileHeader<tensor>(i, fieldName, substep);
        }

        file(i) << endl;
    }
}


Foam::word Foam::functionObjects::residualsBCFDK::subStepName(const label substep)
{
    word substepName = "";
    if (substep>0)
    {
        substepName += "#" + Foam::name(substep);
    }
    return substepName;
}


bool Foam::functionObjects::residualsBCFDK::execute()
{
    return true;
}


bool Foam::functionObjects::residualsBCFDK::write()
{
    logFiles::write();

    if (Pstream::master())
    {
        //Start the line tag
        Log << name() <<  token::COLON << token::SPACE;

        if (writeToFile_)
        {
            writeTime(file(0));
            writeTime(file(1));
        }

        forAll(fieldSet_, fieldi)
        {
            const word& fieldName = fieldSet_[fieldi];
            const label substep = fieldSetSubStep_[fieldi];

            writeResidual<scalar>(fieldName, substep);
            writeResidual<vector>(fieldName, substep);
            writeResidual<sphericalTensor>(fieldName, substep);
            writeResidual<symmTensor>(fieldName, substep);
            writeResidual<tensor>(fieldName, substep);
        }

        Log << nl << endl;

        if (writeToFile_)
        {
            file(0) << endl;
            file(1) << endl;
        }
    }

    return true;
}


// ************************************************************************* //
