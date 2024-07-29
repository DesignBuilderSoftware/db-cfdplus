/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Inherits from OpenFOAM's readFields function object.
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
    Foam::functionObjects::readFieldsBCFDK

Group
    grpFieldFunctionObjects

Description
    Reads fields from the time directories and adds them to the mesh database
    for further post-processing.

    Added \c onConstruction option, because we need it to load right when it
    starts.

    Example of function object specification:
    \verbatim
    readFieldsBCFDK1
    {
        type        readFieldsBCFDK;
        libs        ("libfunctionObjectsBCFDK.so");
        ...
        onConstruction true;
        neverExecute   true;
        fields
        (
            U
            p
        );
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: readFieldsBCFDK   | yes         |
        fields       | list of fields to read  |  no         |
        onConstruction | read when construction is complete | yes |
        neverExecute | do not execute at any time | yes |
    \endtable

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::readFields

SourceFiles
    readFieldsBCFDK.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_readFieldsBCFDK_H
#define functionObjects_readFieldsBCFDK_H

#include "readFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class readFieldsBCFDK Declaration
\*---------------------------------------------------------------------------*/

class readFieldsBCFDK
:
    public readFields
{
protected:

    // Protected data

        //- Whether to read during construction
        Switch onConstruction_;

        //- Whether to read during execution (which will clean up)
        Switch neverExecute_;

public:

    //- Runtime type information
    TypeName("readFieldsBCFDK");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        readFieldsBCFDK
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        readFieldsBCFDK(const readFieldsBCFDK&) = delete;

    //- Destructor
    virtual ~readFieldsBCFDK();


    // Member Functions

        //- Read the additional option
        virtual bool read(const dictionary&);

        //- Read the fields
        virtual bool execute();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const readFieldsBCFDK&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
