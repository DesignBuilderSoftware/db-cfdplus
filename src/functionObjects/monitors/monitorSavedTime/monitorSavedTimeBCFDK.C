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

#include "monitorSavedTimeBCFDK.H"
#include "dictionary.H"
#include "error.H"
#include "OSspecific.H"
#include "addToRunTimeSelectionTable.H"

#if defined( WIN32 ) || defined( WIN64 )
#include "Time.T.H"
#else
#include "Time.H"
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(monitorSavedTimeBCFDK, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        monitorSavedTimeBCFDK,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::monitorSavedTimeBCFDK::monitorSavedTimeBCFDK
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    time_(runTime)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::monitorSavedTimeBCFDK::~monitorSavedTimeBCFDK()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::functionObjects::monitorSavedTimeBCFDK::execute()
{
    // Do nothing - only valid on write
    return true;
}


bool Foam::functionObjects::monitorSavedTimeBCFDK::write()
{
    //Current time
    Info<< "Saved Time = " << time_.timeName()
        //End line twice
        << endl
        << endl;

    return true;
}

// ************************************************************************* //
