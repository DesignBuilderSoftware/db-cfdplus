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
#include "dictionary.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::faceZoneEntry> Foam::faceZoneEntry::New
(
    Istream& is
)
{
    if (debug)
    {
        Info<< "faceZone::New(Istream& is) : constructing faceZone "
            << endl;
    }

    return autoPtr<faceZoneEntry>(new faceZoneEntry(is));
}


// ************************************************************************* //
