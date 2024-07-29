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
    Foam::pointZoneEntry

Description
    A subset of mesh points organised as a primitive patch.

    For quick check whether a point belongs to the zone use the lookup
    mechanism in pointZoneMesh, where all the zoned points are registered
    with their zone number.

SourceFiles
    pointZoneEntry.C
    pointZoneEntryNew.C

\*---------------------------------------------------------------------------*/

#ifndef pointZoneEntry_H
#define pointZoneEntry_H

#include "entry.H"
#include "boolList.H"
#include "primitiveFacePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;

// Forward declaration of friend functions and operators

class pointZoneEntry;
Ostream& operator<<(Ostream&, const pointZoneEntry&);


/*---------------------------------------------------------------------------*\
                          Class pointZoneEntry Declaration
\*---------------------------------------------------------------------------*/

class pointZoneEntry
{

protected:

    // Protected data

        //The loaded-interpreted dictionary
        autoPtr<entry> dictEntry_;

        //- Name of zone
        word name_;

        //- Type of zone
        word type_;

        //- List of indexes of the points assigned to this pointZoneEntry
        labelList pointLabels_;

public:

    // Static data members


    //- Runtime type information
    TypeName("pointZoneEntry");


    // Constructors

        //The only constructor we need, which creates based on the
        //intepretation
        pointZoneEntry(Istream& is);

        pointZoneEntry(const pointZoneEntry&);

        //- Clone
        inline autoPtr<pointZoneEntry> clone() const
        {
          return autoPtr<pointZoneEntry>(new pointZoneEntry(*this));
        }

        //- Construct on freestore from Istream and return
        static autoPtr<pointZoneEntry> New(Istream& is);


    //- Destructor
        virtual ~pointZoneEntry();


    // Member Functions

        //- Write
        virtual void write(Ostream&) const;

        //- Write dictionary
        virtual void writeDict(Ostream&) const;

    // I-O

        //- Ostream Operator
        friend Ostream& operator<<(Ostream&, const pointZoneEntry&);

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const pointZoneEntry&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
