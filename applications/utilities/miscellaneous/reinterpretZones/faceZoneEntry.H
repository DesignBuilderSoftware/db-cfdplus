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
    Foam::faceZoneEntry

Description
    A subset of mesh faces organised as a primitive patch.

    For quick check whether a face belongs to the zone use the lookup
    mechanism in faceZoneMesh, where all the zoned faces are registered
    with their zone number.

SourceFiles
    faceZoneEntry.C
    faceZoneEntryNew.C

\*---------------------------------------------------------------------------*/

#ifndef faceZoneEntry_H
#define faceZoneEntry_H

#include "entry.H"
#include "boolList.H"
#include "primitiveFacePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;

// Forward declaration of friend functions and operators

class faceZoneEntry;
Ostream& operator<<(Ostream&, const faceZoneEntry&);


/*---------------------------------------------------------------------------*\
                          Class faceZoneEntry Declaration
\*---------------------------------------------------------------------------*/

class faceZoneEntry
{

protected:

    // Protected data

        //The loaded-interpreted dictionary
        autoPtr<entry> dictEntry_;

        //- Name of zone
        word name_;

        //- Type of zone
        word type_;

        //- List of indexes of the faces assigned to this faceZoneEntry
        labelList faceLabels_;

        //- Flip map for all faces in the zone.  Set to true if the
        //  face needs to be flipped to achieve the correct orientation.
        boolList flipMap_;

public:

    // Static data members


    //- Runtime type information
    TypeName("faceZoneEntry");


    // Constructors

        //The only constructor we need, which creates based on the
        //intepretation
        faceZoneEntry(Istream& is);

        faceZoneEntry(const faceZoneEntry&);

        //- Clone
        inline autoPtr<faceZoneEntry> clone() const
        {
          return autoPtr<faceZoneEntry>(new faceZoneEntry(*this));
        }

        //- Construct on freestore from Istream and return
        static autoPtr<faceZoneEntry> New(Istream& is);


    //- Destructor
        virtual ~faceZoneEntry();


    // Member Functions

        //- Write
        virtual void write(Ostream&) const;

        //- Write dictionary
        virtual void writeDict(Ostream&) const;

    // I-O

        //- Ostream Operator
        friend Ostream& operator<<(Ostream&, const faceZoneEntry&);

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const faceZoneEntry&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
