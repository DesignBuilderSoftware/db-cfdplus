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
    Foam::faceEntryList

Description
    A list of face zones.

SourceFiles
    faceEntryList.C

\*---------------------------------------------------------------------------*/

#ifndef faceEntryList_H
#define faceEntryList_H

#if defined( WIN32 ) || defined( WIN64 )
#include "PtrList.T.H"
#else
#include "PtrList.H"
#endif

#include "faceZoneEntry.H"
//#include "entry.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class entryList Declaration
\*---------------------------------------------------------------------------*/

class faceEntryList
:
    public regIOobject,
    public PtrList<faceZoneEntry>
{

public:

    // Constructors

        //- Read constructor given IOobject and a MeshType reference
        faceEntryList
        (
            const IOobject&
        );

        //- Construct given size
        faceEntryList
        (
            const IOobject&,
            const label size
        );

        //- Disallow construct as copy
        faceEntryList(const faceEntryList&) = delete;


    //- Destructor
    ~faceEntryList();


    // Member Functions

        //- writeData member function required by regIOobject
        bool writeData(Ostream&) const;

    // Member Operators

        //- Disallow assignment
        void operator=(const faceEntryList&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
