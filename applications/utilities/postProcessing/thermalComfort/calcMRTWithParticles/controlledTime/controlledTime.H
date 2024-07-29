/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Wrapper class for OpenFOAM's Time class.
    Developed for blueCFD(R)-Kernel
    Copyright (C) 2016 FSD blueCAPE Lda  http://www.bluecape.com.pt
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
    Foam::controlledTime

Description
    Class that derives from Foam::Time and provides more control over time
    parameters, such as "writeInternal".

SourceFiles
    controlledTime.C

\*---------------------------------------------------------------------------*/

#ifndef controlledTime_H
#define controlledTime_H

#if defined(WIN32) || defined(WIN64)
#include "Time.T.H"
#else
#include "Time.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class controlledTime Declaration
\*---------------------------------------------------------------------------*/

class controlledTime
:
    public Time
{
public:

    TypeName("controlledTime");

    // Constructors

        //- Construct given name of dictionary to read and argument list
        controlledTime
        (
            const word& name,
            const argList& args,
            const word& systemName = "system",
            const word& constantName = "constant"
        );

        //- Construct given name of dictionary to read, rootPath and casePath
        controlledTime
        (
            const word& name,
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName = "system",
            const word& constantName = "constant",
            const bool enableFunctionObjects = true
        );

        //- Construct given dictionary, rootPath and casePath
        controlledTime
        (
            const dictionary& dict,
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName = "system",
            const word& constantName = "constant",
            const bool enableFunctionObjects = true
        );

        //- Construct given endTime, rootPath and casePath
        controlledTime
        (
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName = "system",
            const word& constantName = "constant",
            const bool enableFunctionObjects = true
        );


    //- Destructor
    virtual ~controlledTime();


    // Member Functions

        // Access

            //- Return true if time currently being sub-cycled, otherwise false
            scalar writeInterval() const
            {
                return writeInterval_;
            }

        // Edit

            //- Reset the time and time-index to those of the given time
            virtual void setWriteInterval(const scalar);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
