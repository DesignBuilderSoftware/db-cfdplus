/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on blueCFD-Kernel's writeNowAtEnd function object.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2013-2018 FSD blueCAPE Lda  http://www.bluecape.com.pt
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
    Foam::functionObjects::monitorSavedTimeBCFDK

Group
    grpUtilitiesFunctionObjects

Description
    This function object outputs the current time.

    Example of function object specification:
    \verbatim
    savedTime
    {
        type        monitorSavedTimeBCFDK;
        functionObjectLibs ("libMonitorsBCFDK.so");
        ...
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: monitorSavedTimeBCFDK  | yes         |
    \endtable

See also
    Foam::functionObject
    Foam::OutputFilterFunctionObject

SourceFiles
    monitorSavedTimeBCFDK.C

\*---------------------------------------------------------------------------*/

#ifndef monitorSavedTimeBCFDK_H
#define monitorSavedTimeBCFDK_H

#include "functionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class monitorSavedTimeBCFDK Declaration
\*---------------------------------------------------------------------------*/

class monitorSavedTimeBCFDK
:
    public functionObject
{
    // Private Data

        //- Reference to the Time
        const Time& time_;

public:

    //- Runtime type information
    TypeName("monitorSavedTimeBCFDK");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        monitorSavedTimeBCFDK
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        monitorSavedTimeBCFDK(const monitorSavedTimeBCFDK&) = delete;

    //- Destructor
    virtual ~monitorSavedTimeBCFDK();


    // Member Functions

        //- Do nothing
        virtual bool execute();

        //- Write the current time
        virtual bool write();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const monitorSavedTimeBCFDK&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
