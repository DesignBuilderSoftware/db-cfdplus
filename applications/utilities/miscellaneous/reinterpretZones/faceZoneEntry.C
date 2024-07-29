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

#include "faceZoneEntry.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "demandDrivenData.H"
#include "mapPolyMesh.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceZoneEntry, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceZoneEntry::faceZoneEntry(Istream& is) :
  dictEntry_(entry::New(is)),
  name_(dictEntry_().keyword()),
  type_(dictEntry_().dict().lookup("type")),
  faceLabels_(dictEntry_().dict().lookup("faceLabels")),
  flipMap_(dictEntry_().dict().lookup("flipMap"))
{
}


Foam::faceZoneEntry::faceZoneEntry(const faceZoneEntry &other) :
  dictEntry_(), //nothing to do on this one
  name_(other.name_),
  type_(other.type_),
  faceLabels_(other.faceLabels_),
  flipMap_(other.flipMap_)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faceZoneEntry::~faceZoneEntry()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faceZoneEntry::write(Ostream& os) const
{
    writeDict(os);
}


void Foam::faceZoneEntry::writeDict(Ostream& os) const
{
    os  << nl << name_ << nl << token::BEGIN_BLOCK << nl
        << "    type " << type_ << token::END_STATEMENT << nl;

    writeEntry(os, "faceLabels", faceLabels_);
    writeEntry(os, "flipMap", flipMap_);

    os  << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const faceZoneEntry& zn)
{
    zn.write(os);
    os.check("Ostream& operator<<(Ostream&, const faceZoneEntry&");
    return os;
}


// ************************************************************************* //
