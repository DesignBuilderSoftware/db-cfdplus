/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Derived from OpenFOAM's volFieldValue function object.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2013-2019 FSD blueCAPE Lda  http://www.bluecape.com.pt
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
    Foam::functionObjects::fieldValues::volumeFieldValueBCFDK

Group
    grpFieldFunctionObjects

Description
    Provides a 'volRegion' specialization of the fieldValue function object.

    Given a list of user-specified fields and a 'volRegion', a number of
    operations can be performed, such as sums, averages and integrations.

    Examples of function object specification:
    \verbatim
    volFieldValue1
    {
        type            volFieldValueBCFDK;
        libs            ("libMonitorsBCFDK.so");

        log             true;
        writeControl    writeTime;
        writeFields     true;

        regionType      cellZone;
        name            c0;
        operation       volAverage;

        weightField     alpha1;

        fields
        (
            p
            U
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

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: volFieldValue            | yes   |
        log          | Write data to standard output       | no    | no
        writeFields  | Write the region field values       | yes   |
        regionType   | volRegion type: see below           | yes   |
        name         | Name of volRegion if required       | no    |
        operation    | Operation to perform                | yes   |
        weightField  | Name of field to apply weighting    | no    |
        fields       | List of fields to operate on        | yes   |
        fields       | list of fields to operate on | yes    |
        identifiers  | corresponding list of names for each field entry | yes |
        modes        | list of calculation modes:
                       scalar, magnitude, componentX,
                       componentY, componentZ | yes | magnitude
    \endtable

    Where \c regionType is defined by
    \plaintable
        cellZone     | requires a 'name' entry to specify the cellZone
        all          | all cells
    \endplaintable

    The \c operation is one of:
    \plaintable
       none                 | No operation
       sum                  | Sum
       weightedSum          | Weighted sum
       sumMag               | Sum of component magnitudes
       average              | Ensemble average
       weightedAverage      | Weighted average
       volAverage           | Volume weighted average
       weightedVolAverage   | Weighted volume average
       volIntegrate         | Volume integral
       weightedVolIntegrate | Weighted volume integral
       min                  | Minimum
       max                  | Maximum
       CoV                  | Coefficient of variation: standard deviation/mean
    \endplaintable

See also
    Foam::functionObjects::fieldValues::fieldValue
    Foam::functionObjects::volRegion
    Foam::functionObject

SourceFiles
    volumeFieldValueBCFDK.C
    volumeFieldValueBCFDKTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_volumeFieldValueBCFDK_H
#define functionObjects_volumeFieldValueBCFDK_H

#include "volFieldValue.H"
#include "commonDictBCFDK.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace fieldValues
{

/*---------------------------------------------------------------------------*\
                      Class volumeFieldValueBCFDK Declaration
\*---------------------------------------------------------------------------*/

class volumeFieldValueBCFDK
:
    public volFieldValue,
    public commonDictBCFDK
{

public:

    //- Run-time type information
    TypeName("volFieldValueBCFDK");


    // Constructors

        //- Construct from name, Time and dictionary
        volumeFieldValueBCFDK
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct from name, objectRegistry and dictionary
        volumeFieldValueBCFDK
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );


    //- Destructor
    virtual ~volumeFieldValueBCFDK();


    // Public Member Functions

        //- Custom templated helper function to output field values
        template<class Type>
        bool writeValuesBCFDK
        (
            const word& fieldName,
            const modeType& mode,
            const word& identifier
        );

        //- Read from dictionary
        virtual bool read(const dictionary&);

        //- Calculate and write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fieldValues
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "volumeFieldValueBCFDKTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
