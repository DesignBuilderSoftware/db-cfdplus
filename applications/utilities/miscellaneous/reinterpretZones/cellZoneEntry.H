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

Class
    Foam::cellZoneEntry

Description
    A subset of mesh cells organised as a primitive patch.

    For quick check whether a cell belongs to the zone use the lookup
    mechanism in cellZoneMesh, where all the zoned cells are registered
    with their zone number.

SourceFiles
    cellZoneEntry.C
    cellZoneEntryNew.C

\*---------------------------------------------------------------------------*/

#ifndef cellZoneEntry_H
#define cellZoneEntry_H

#include "entry.H"
#include "boolList.H"
#include "primitiveFacePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;

// Forward declaration of friend functions and operators

class cellZoneEntry;
Ostream& operator<<(Ostream&, const cellZoneEntry&);


/*---------------------------------------------------------------------------*\
                          Class cellZoneEntry Declaration
\*---------------------------------------------------------------------------*/

class cellZoneEntry
{

protected:

    // Protected data

        //The loaded-interpreted dictionary
        autoPtr<entry> dictEntry_;

        //- Name of zone
        word name_;

        //- Type of zone
        word type_;

        //- List of indexes of the cells assigned to this cellZoneEntry
        labelList cellLabels_;


public:

    // Static data members


    //- Runtime type information
    TypeName("cellZoneEntry");


    // Constructors

        //The only constructor we need, which creates based on the
        //intepretation
        cellZoneEntry(Istream& is);

        cellZoneEntry(const cellZoneEntry&);

        //- Clone
        inline autoPtr<cellZoneEntry> clone() const
        {
          return autoPtr<cellZoneEntry>(new cellZoneEntry(*this));
        }

        //- Construct on freestore from Istream and return
        static autoPtr<cellZoneEntry> New(Istream& is);


    //- Destructor
        virtual ~cellZoneEntry();


    // Member Functions

        //- Write
        virtual void write(Ostream&) const;

        //- Write dictionary
        virtual void writeDict(Ostream&) const;

    // I-O

        //- Ostream Operator
        friend Ostream& operator<<(Ostream&, const cellZoneEntry&);

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cellZoneEntry&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
