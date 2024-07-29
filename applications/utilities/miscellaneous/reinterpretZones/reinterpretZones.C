/*---------------------------------------------------------------------------*\
Changes, Authors and Copyright
    Based on various OpenFOAM utilities.
    Developed for blueCFD(R)-Kernel
    Copyright (C) 2014-2018 FSD blueCAPE Lda  http://www.bluecape.com.pt
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

Application
    reinterpretZones

Description
    Read the zones dictionaries and write the dictionary in ascii.

Usage
    - reinterpretFaceZones inputDictFolder [OPTION]

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "faceEntryList.H"
#include "cellEntryList.H"
#include "pointEntryList.H"
#include "polyMesh.H"

#if defined( WIN32 ) || defined( WIN64 )
#include "Time.T.H"
#else
#include "Time.H"
#endif

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("inputDictFolder");
    #include "setRootCase.H"
    #include "createTime.H"

    const string dictFolderName = args[1];

    IOobject faceZonesHeader
    (
        dictFolderName + "/faceZones",
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if(faceZonesHeader.typeHeaderOk<regIOobject>(true))
    {
        Info<< "Processing " << faceZonesHeader.name() << endl;

        faceEntryList faceDict(faceZonesHeader);

        faceDict.writeObject(
            IOstream::ASCII,
            IOstream::currentVersion,
            IOstream::UNCOMPRESSED,
            true
            );
    }
    else
    {
        Info<< "Skipping " << faceZonesHeader.name() << endl;
    }


    IOobject cellZonesHeader
    (
        dictFolderName + "/cellZones",
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if(cellZonesHeader.typeHeaderOk<regIOobject>(true))
    {
        Info<< "Processing " << cellZonesHeader.name() << endl;

        cellEntryList cellDict(cellZonesHeader);

        cellDict.writeObject(
            IOstream::ASCII,
            IOstream::currentVersion,
            IOstream::UNCOMPRESSED,
            true
            );
    }
    else
    {
        Info<< "Skipping " << cellZonesHeader.name() << endl;
    }


    IOobject pointZonesHeader
    (
        dictFolderName + "/pointZones",
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if(pointZonesHeader.typeHeaderOk<regIOobject>(true))
    {
        Info<< "Processing " << pointZonesHeader.name() << endl;

        pointEntryList pointDict(pointZonesHeader);

        pointDict.writeObject(
            IOstream::ASCII,
            IOstream::currentVersion,
            IOstream::UNCOMPRESSED,
            true
            );
    }
    else
    {
        Info<< "Skipping " << pointZonesHeader.name() << endl;
    }

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
