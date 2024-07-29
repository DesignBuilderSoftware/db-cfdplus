/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's fieldMinMax function object.
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
    Foam::functionObjects::fieldMinMaxBCFDK

Group
    grpFieldFunctionObjects

Description
    Calculates the value of scalar minimum and maximum for a list of
    user-specified fields.

    For variables with a rank greater than zero, either the min/max of a
    component value or the magnitude is reported.  When operating in parallel,
    the processor owning the value is also given.

    The additional features are regarding the \c dataMode and \c identifiers,
    because more fine grained information is needed for blueCFD-Kernel. Also
    removed the report of the location of said value.

    Example of function object specification:
    \verbatim
    fieldMinMaxBCFDK1
    {
        type        fieldMinMaxBCFDK;
        functionObjectLibs ("libMonitorsBCFDK.so");
        ...
        operation   maximum;
        dataMode    cells;
        fields
        (
            U
            p
        );
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
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: fieldMinMaxBCFDK  | yes         |
        operation    | calculation modes: maximum or minimum | yes | maximum
        dataMode     | retrieve data from: cells, boundaries or both
                                               | yes | both
        identifiers  | corresponding list of names for each field entry | yes |
        modes        | list of calculation modes: scalar, magnitude,
                       componentX, componentY, componentZ | yes | magnitude
    \endtable

    Output data is written to the file \<timeDir\>/fieldMinMax.dat

Note
    - Based on OpenFOAM's fieldMinMax function object.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    fieldMinMaxBCFDK.C

\*---------------------------------------------------------------------------*/

#ifndef fieldMinMaxBCFDK_H
#define fieldMinMaxBCFDK_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "vector.H"
#include "commonDictBCFDK.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class fieldMinMaxBCFDK Declaration
\*---------------------------------------------------------------------------*/

class fieldMinMaxBCFDK
:
    public fvMeshFunctionObject,
    public logFiles,
    public commonDictBCFDK
{
public:

    enum operationType
    {
        mdMaximum,
        mdMinimum
    };

    enum dataModeType
    {
        mdCells,
        mdBoundaries,
        mdBoth
    };

protected:

    // Protected data

        //- Operation type names
        static const NamedEnum<operationType, 2> operationTypeNames_;

        //- Data operation type names
        static const NamedEnum<dataModeType, 3> dataModeTypeNames_;

        //- Operation mode (minimum or maximum)
        operationType operation_;

        //- Operation mode (cells, boundaries or both)
        dataModeType dataMode_;

        //- Fields to assess min/max
        wordList fieldSet_;

    // Protected Member Functions

        //- Calculate the field min/max
        template<class Type>
        void calcMinMaxFields
        (
            const word& fieldName,
            const modeType& mode,
            const word& append2FieldName
        );

        //- Write file header, implementing from logFiles
        void writeFileHeader(const label i);

public:

    //- Runtime type information
    TypeName("fieldMinMaxBCFDK");


    // Constructors

        //- Construct from Time and dictionary
        fieldMinMaxBCFDK
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        fieldMinMaxBCFDK(const fieldMinMaxBCFDK&) = delete;


    //- Destructor
    virtual ~fieldMinMaxBCFDK();


    // Member Functions

        //- Read the field min/max data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the fieldMinMaxBCFDK
        virtual bool write();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fieldMinMaxBCFDK&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fieldMinMaxBCFDKTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
