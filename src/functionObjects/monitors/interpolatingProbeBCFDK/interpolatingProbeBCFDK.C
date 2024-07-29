/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's surfaceFieldValue, fieldMinMax and probe function
    objects.
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

#include "interpolatingProbeBCFDK.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(interpolatingProbesBCFDK, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        interpolatingProbesBCFDK,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::interpolatingProbesBCFDK::findElement
(
    const fvMesh& mesh
)
{
    if (debug)
    {
        Info<< "resetting sample location" << endl;
    }

    // Reset these values just in case
    elementProbe_ = -1;
    faceProbe_ = -1;

    label celli = mesh.findCell(probePosition_);
    elementProbe_ = celli;

    if (celli != -1)
    {
        const labelList& cellFaces = mesh.cells()[celli];
        const vector& cellCentre = mesh.cellCentres()[celli];
        scalar minDistance = GREAT;
        label minFaceID = -1;
        forAll(cellFaces, i)
        {
            label facei = cellFaces[i];
            vector dist = mesh.faceCentres()[facei] - cellCentre;
            if (mag(dist) < minDistance)
            {
                minDistance = mag(dist);
                minFaceID = facei;
            }
        }
        faceProbe_ = minFaceID;
    }
    else
    {
        faceProbe_ = -1;
    }

    if (debug && (elementProbe_ != -1 || faceProbe_ != -1))
    {
        Pout<< "found point " << probePosition_
            << " in cell " << elementProbe_
            << " and face " << faceProbe_ << endl;
    }

    //label celli = elementProbe_;
    label facei = faceProbe_;

    // Check at least one processor with cell.
    reduce(celli, maxOp<label>());
    reduce(facei, maxOp<label>());

    if (celli == -1)
    {
        if (Pstream::master())
        {
            WarningInFunction
                << "Did not find location " << probePosition_
                << " in any cell. Skipping location." << endl;
        }
    }
    else if (facei == -1)
    {
        if (Pstream::master())
        {
            WarningInFunction
                << "Did not find location " << probePosition_
                << " in any face. Skipping location." << endl;
        }
    }
    else
    {
        // Make sure location not on two domains.
        if (elementProbe_ != -1 && elementProbe_ != celli)
        {
            WarningInFunction
                << "Location " << probePosition_
                << " seems to be on multiple domains:"
                << " cell " << elementProbe_
                << " on my domain " << Pstream::myProcNo()
                    << " and cell " << celli << " on some other domain."
                << endl
                << "This might happen if the probe location is on"
                << " a processor patch. Change the location slightly"
                << " to prevent this." << endl;
        }

        if (faceProbe_ != -1 && faceProbe_ != facei)
        {
            WarningInFunction
                << "Location " << probePosition_
                << " seems to be on multiple domains:"
                << " cell " << faceProbe_
                << " on my domain " << Pstream::myProcNo()
                    << " and face " << facei << " on some other domain."
                << endl
                << "This might happen if the probe location is on"
                << " a processor patch. Change the location slightly"
                << " to prevent this." << endl;
        }
    }
}

void Foam::functionObjects::interpolatingProbesBCFDK::writeFileHeader
(
    const label i
)
{
    commonDictBCFDK::writeFileHeader(file(i), "interpolatingProbe");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::interpolatingProbesBCFDK::
interpolatingProbesBCFDK
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    commonDictBCFDK(dict)
{
    logFiles::resetName(typeName);

    read(dict);
}

Foam::functionObjects::interpolatingProbesBCFDK::
interpolatingProbesBCFDK
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, obr, dict),
    logFiles(obr_, name),
    commonDictBCFDK(dict)
{
    logFiles::resetName(typeName);

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::interpolatingProbesBCFDK::
~interpolatingProbesBCFDK()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::interpolatingProbesBCFDK::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);
    commonDictBCFDK::read(dict);

    debug = dict.lookupOrDefault<label>("debug", debug);

    dict.lookup("fields") >> fields_;
    dict.lookup("probePosition") >> probePosition_;
    dict.lookup("interpolationScheme") >> interpolationScheme_;

    // Initialise cell to sample from supplied location
    findElement(mesh_);

    return true;
}


bool Foam::functionObjects::interpolatingProbesBCFDK::execute()
{
    // Do nothing here
    return true;
}


bool Foam::functionObjects::interpolatingProbesBCFDK::write()
{
    //Start the line tag
    Log << name() <<  token::COLON << token::SPACE;

    if (logToFile_ && Pstream::master())
    {
        logFiles::write();
        writeTime(file());
    }

    // Process the fields
    forAll(fields_, i)
    {
        const word& fieldName = fields_[i];
        const modeType& modeName = modes_[i];
        const word& idName = identifiers_[i];
        bool ok = false;

        ok = ok || probePosition<scalar>
            (
                fieldName,
                modeName,
                idName
            );
        ok = ok || probePosition<vector>
            (
                fieldName,
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

    //End line twice
    Log << nl << endl;

    return true;
}


// ************************************************************************* //
