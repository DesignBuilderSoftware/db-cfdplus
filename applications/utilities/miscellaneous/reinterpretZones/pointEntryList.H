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
    Foam::pointEntryList

Description
    A list of point zones.

SourceFiles
    pointEntryList.C

\*---------------------------------------------------------------------------*/

#ifndef pointEntryList_H
#define pointEntryList_H

#if defined( WIN32 ) || defined( WIN64 )
#include "PtrList.T.H"
#else
#include "PtrList.H"
#endif

#include "pointZoneEntry.H"
//#include "entry.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pointList Declaration
\*---------------------------------------------------------------------------*/

class pointEntryList
:
    public regIOobject,
    public PtrList<pointZoneEntry>
{

public:

    // Constructors

        //- Read constructor given IOobject and a MeshType reference
        pointEntryList
        (
            const IOobject&
        );

        //- Construct given size
        pointEntryList
        (
            const IOobject&,
            const label size
        );

        //- Disallow construct as copy
        pointEntryList(const pointEntryList&) = delete;

    //- Destructor
    ~pointEntryList();


    // Member Functions

        //- writeData member function required by regIOobject
        bool writeData(Ostream&) const;

    // Member Operators

        //- Disallow assignment
        void operator=(const pointEntryList&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
