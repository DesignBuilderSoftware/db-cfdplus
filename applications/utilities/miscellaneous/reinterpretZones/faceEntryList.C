/*---------------------------------------------------------------------------*\
Changes, Authors and Copyright
    Developed for blueCFD(R)-Kernel
    Copyright (C) 2014-2016 FSD blueCAPE Lda  http://www.bluecape.com.pt
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

#include "faceEntryList.H"
#include "demandDrivenData.H"
#include "stringListOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Read constructor given IOobject and a MeshType reference
Foam::faceEntryList::faceEntryList
(
    const IOobject& io
)
:
    regIOobject(io),
    PtrList<faceZoneEntry>(readStream(typeName))
{
    close();
}


// Construct given size. Zones will be set later
Foam::faceEntryList::faceEntryList
(
    const IOobject& io,
    const label size
)
:
    regIOobject(io),
    PtrList<faceZoneEntry>(size)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faceEntryList::~faceEntryList()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// writeData member function required by regIOobject
bool Foam::faceEntryList::writeData(Ostream& os) const
{
    os  << *this;
    return os.good();
}

// ************************************************************************* //
