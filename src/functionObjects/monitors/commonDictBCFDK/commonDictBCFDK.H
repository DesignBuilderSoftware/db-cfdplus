/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Derived from OpenFOAM's class structures.
    Developed for blueCFD(R)-Kernel
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
    Foam::commonDictBCFDK

Group
    grpFieldFunctionObjects

Description
    Provides a

    Example of the additional dictionary entries that have been added:
    \verbatim
    <block name here>
    {
        //...
        identifiers
        (
            "U|Mag"
            "Pressure"
        );
        modes
        (
            scalar
            magnitude
        );
        //...
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        identifiers  | corresponding list of names for each field entry | yes |
        modes        | list of calculation modes:
                       scalar, magnitude, componentX,
                       componentY, componentZ | yes | magnitude
    \endtable

SourceFiles
    commonDictBCFDK.C

\*---------------------------------------------------------------------------*/

#ifndef commonDictBCFDK_H
#define commonDictBCFDK_H

#if defined( WIN32 ) || defined( WIN64 )
#include "NamedEnum.T.H"
#include "Field.T.H"
#else
#include "NamedEnum.H"
#include "Field.H"
#endif

#include "OFstream.H"
#include "dictionary.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class commonDictBCFDK Declaration
\*---------------------------------------------------------------------------*/

class commonDictBCFDK
{

public:

    // Public data types

        //- Mode type enumeration
        enum modeType
        {
            mdScalar,
            mdMag,
            mdCmptX,
            mdCmptY,
            mdCmptZ,
        };

        //- Mode type names
        static const NamedEnum<modeType, 5> modeTypeNames_;

    // Protected Member Functions

        template<class ReturnType, class Type>
        tmp<ReturnType> convertField
        (
            const Type& field,
            const modeType& mode
        ) const;

protected:

    // Protected data

        //- List of text strings to be appended to the name on the report
        wordList identifiers_;

        //- Modes for which component to get
        List<modeType> modes_;

        //- Log to file as well
        Switch logToFile_;

public:

    // Constructors

        //- Construct from dictionary
        commonDictBCFDK
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~commonDictBCFDK();


    // Public Member Functions

        //- Read from dictionary
        virtual bool read(const dictionary&);

        //- Write common header to file
        void writeFileHeader
        (
            OFstream& fileOFS,
            const word commonOperation
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "commonDictBCFDKTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
