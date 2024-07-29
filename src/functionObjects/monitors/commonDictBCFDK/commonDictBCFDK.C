/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Derived from OpenFOAM's class structures.
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

#include "commonDictBCFDK.H"
#include "wordList.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* Foam::NamedEnum
<
    Foam::commonDictBCFDK::modeType,
    5
>::names[] =
{
    "scalar",
    "magnitude",
    "componentX",
    "componentY",
    "componentZ"
};

const Foam::NamedEnum
<
    Foam::commonDictBCFDK::modeType,
    5
>
Foam::commonDictBCFDK::modeTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::commonDictBCFDK::
commonDictBCFDK
(
    const dictionary& dict
)
:
    identifiers_(),
    modes_(),
    logToFile_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::commonDictBCFDK::
~commonDictBCFDK()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::commonDictBCFDK::read
(
    const dictionary& dict
)
{
    dict.lookup("identifiers") >> identifiers_;

    wordList modesAsList;
    dict.lookup("modes") >> modesAsList;

    modes_.setSize(modesAsList.size(), mdScalar);
    forAll(modesAsList, modeI)
    {
        modes_[modeI] = modeTypeNames_[modesAsList[modeI]];
    }

    logToFile_ = dict.lookupOrDefault("logToFile", false);

    return true;
}

void Foam::commonDictBCFDK::writeFileHeader
(
    OFstream& fileOFS,
    const word commonOperation
)
{
    fileOFS << "# Time";

    forAll(identifiers_, idi)
    {
        fileOFS
            << tab
            << commonOperation
            << "(" << identifiers_[idi] << ")";
    }

    fileOFS << endl;
}
// ************************************************************************* //
