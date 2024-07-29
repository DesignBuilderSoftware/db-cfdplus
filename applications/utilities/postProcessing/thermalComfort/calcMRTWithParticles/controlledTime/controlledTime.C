/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Wrapper class for OpenFOAM's controlledTime class.
    Developed for blueCFD(R)-Kernel
    Copyright (C) 2016 FSD blueCAPE Lda  http://www.bluecape.com.pt
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

#include "controlledTime.H"
#include "argList.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(controlledTime, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::controlledTime::controlledTime
(
    const word& controlDictName,
    const fileName& rootPath,
    const fileName& caseName,
    const word& systemName,
    const word& constantName,
    const bool enableFunctionObjects
)
:
    Time
    (
        controlDictName,
        rootPath,
        caseName,
        systemName,
        constantName,
        enableFunctionObjects
    )
{
}


Foam::controlledTime::controlledTime
(
    const word& controlDictName,
    const argList& args,
    const word& systemName,
    const word& constantName
)
:
    Time
    (
        controlDictName,
        args,
        systemName,
        constantName
    )
{
}


Foam::controlledTime::controlledTime
(
    const dictionary& dict,
    const fileName& rootPath,
    const fileName& caseName,
    const word& systemName,
    const word& constantName,
    const bool enableFunctionObjects
)
:
    Time
    (
        dict,
        rootPath,
        caseName,
        systemName,
        constantName,
        enableFunctionObjects
    )
{
}


Foam::controlledTime::controlledTime
(
    const fileName& rootPath,
    const fileName& caseName,
    const word& systemName,
    const word& constantName,
    const bool enableFunctionObjects
)
:
    Time
    (
        rootPath,
        caseName,
        systemName,
        constantName,
        enableFunctionObjects
    )
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::controlledTime::~controlledTime()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::controlledTime::setWriteInterval
(
    const scalar newWriteInterval
)
{
    writeInterval_ = newWriteInterval;
}


// ************************************************************************* //
