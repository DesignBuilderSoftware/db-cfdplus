/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Inherits from OpenFOAM's readFields function object.
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

#include "readFieldsBCFDK.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(readFieldsBCFDK, 0);
    addToRunTimeSelectionTable(functionObject, readFieldsBCFDK, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::readFieldsBCFDK::readFieldsBCFDK
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    readFields(name, runTime, dict)
{
    read(dict);

    if(onConstruction_)
    {
        readFields::execute();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::readFieldsBCFDK::~readFieldsBCFDK()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::readFieldsBCFDK::read(const dictionary& dict)
{
    readFields::read(dict);

    dict.lookup("onConstruction") >> onConstruction_;
    dict.lookup("neverExecute") >> neverExecute_;

    return true;
}

bool Foam::functionObjects::readFieldsBCFDK::execute()
{
    if(!neverExecute_)
    {
        readFields::execute();
    }

    return true;
}

// ************************************************************************* //
