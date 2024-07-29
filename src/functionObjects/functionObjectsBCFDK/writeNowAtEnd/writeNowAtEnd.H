/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Derives from OpenFOAM's abort function object.
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

Class
    Foam::functionObjects::writeNowAtEnd

Group
    grpUtilitiesFunctionObjects

Description
    Watches for presence of the named file in the $FOAM_CASE directory
    and writeNowAtEnds the calculation if it is present.

    Currently the following action types are supported:
    - noWriteNow
    - writeNow
    - nextWrite

SourceFiles
    writeNowAtEnd.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_writeNowAtEnd_H
#define functionObjects_writeNowAtEnd_H

#include "functionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                      Class writeNowAtEnd Declaration
\*---------------------------------------------------------------------------*/

class writeNowAtEnd
:
    public functionObject
{
    // Private Data

        //- Reference to the Time
        const Time& time_;


public:

    //- Runtime type information
    TypeName("writeNowAtEnd");


    // Constructors

        //- Construct from Time and dictionary
        writeNowAtEnd
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        writeNowAtEnd(const writeNowAtEnd&) = delete;

    //- Destructor
    virtual ~writeNowAtEnd();


    // Member Functions

        //- Do nothing
        virtual bool execute();

        //- Do nothing
        virtual bool write();

        //- Execute at the final time-loop, used for cleanup
        virtual bool end();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const writeNowAtEnd&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
