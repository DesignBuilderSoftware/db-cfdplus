/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Derives from OpenFOAM's abort function object.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2018 FSD blueCAPE Lda  http://www.bluecape.com.pt
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

#include "writeNowAtEnd.H"
#include "dictionary.H"
#include "error.H"
#include "OSspecific.H"
#include "addToRunTimeSelectionTable.H"

#if defined( WIN32 ) || defined( WIN64 )
#include "Time.T.H"
#include "PstreamReduceOps.T.H"
#else
#include "Time.H"
#include "PstreamReduceOps.H"
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(writeNowAtEnd, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        writeNowAtEnd,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::writeNowAtEnd::writeNowAtEnd
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

Foam::functionObjects::writeNowAtEnd::~writeNowAtEnd()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::writeNowAtEnd::execute()
{
    return true;
}


bool Foam::functionObjects::writeNowAtEnd::write()
{
    return true;
}


bool Foam::functionObjects::writeNowAtEnd::end()
{
    if (!time_.writeTime())
    {
        Info<< type() << " has now forced 'writeNow' policy (timeIndex="
            << time_.timeIndex()
            << "): stop+write data"
            << endl;

        //Need non-const permissions, so we will hack our way in
        Time& overrideTime = const_cast<Time&>(time_);

        overrideTime.writeNow();
        overrideTime.functionObjects().execute();
    }

    return true;
}

// ************************************************************************* //
