/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Derived from OpenFOAM's surfaceFieldValue function object.
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
    Foam::functionObjects::fieldValues::surfaceFieldValueBCFDK

Description
    Provides a 'face regionType' variant of the fieldValues function object.

    Given a list of user-specified fields and a selection of mesh (or general
    surface) faces, a number of operations can be performed, such as sums,
    averages and integrations.

    For example, to calculate the volumetric or mass flux across a patch,
    apply the 'sum' operator to the flux field (typically \c phi).

    Examples of function object specification:
    \verbatim
    surfaceFieldValue1
    {
        type            surfaceFieldValueBCFDK;
        libs            ("libMonitorsBCFDK.so");

        log             true; //Should always be set to 'true'
        writeControl    writeTime;
        writeFields     false; //Do not use any 'true'

        surfaceFormat   none; //Do not use any other
        regionType      patch; //faceZones are not interpolated
        name            f0;

        operation       sum;

        weightField     phi;

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
        type         | type name: surfaceFieldValueBCFDK   | yes         |
        log          | write data to standard output | no    | no
        writeFields  | Write the region field values  | yes     |
        writeArea    | Write the area of the selected surface | no |
        surfaceFormat | output value format    | no          |
        regionType   | face regionType: see below  | yes         |
        name         | name of face regionType if required  | no |
        operation    | operation to perform    | yes         |
        weightField  | name of field to apply weighting | no |
        orientedWeightField  | name of oriented field to apply weighting | no |
        scaleFactor  | scale factor            | no          | 1
        fields       | list of fields to operate on | yes    |
        orientedFields | list of oriented fields to operate on | no |
        identifiers  | corresponding list of names for each field entry | yes |
        modes        | list of calculation modes:
                       scalar, magnitude, componentX,
                       componentY, componentZ | yes | magnitude
    \endtable

    Where \c regionType is defined by
    \plaintable
        faceZone     | requires a 'name' entry to specify the faceZone
        patch        | requires a 'name' entry to specify the patch
        sampledSurface | requires a 'sampledSurfaceDict' sub-dictionary
    \endplaintable

    The \c operation is one of:
    \plaintable
       none          | no operation
       sum           | sum
       weightedSum   | weighted sum
       sumMag        | sum of component magnitudes
       sumDirection  | sum values which are positive in given direction
       sumDirectionBalance | sum of balance of values in given direction
       average       | ensemble average
       weightedAverage | weighted average
       areaAverage   | area weighted average
       weightedAreaAverage | weighted area average
       areaIntegrate | area integral
       weightedAreaIntegrate | weighted area integral
       min           | minimum
       max           | maximum
       CoV           | coefficient of variation: standard deviation/mean
       areaNormalAverage| area weighted average in face normal direction
       areaNormalIntegrate | area weighted integral in face normal directon
    \endplaintable

Note
    - The values reported by the areaNormalAverage and areaNormalIntegrate
      operations are written as the first component of a field with the same
      rank as the input field.
    - faces on empty patches get ignored
    - if the field is a volField the \c faceZone can only consist of boundary
      faces
    - the `oriented' entries relate to mesh-oriented fields, such as the
      flux, phi.  These fields will be oriented according to the face normals.
    - using \c sampledSurface:
        - not available for surface fields
        - if interpolate=true they use \c interpolationCellPoint
          otherwise they use cell values
        - each triangle in \c sampledSurface is logically only in one cell
          so interpolation will be wrong when triangles are larger than
          cells.  This can only happen for sampling on a \c triSurfaceMesh
        - take care when using isoSurfaces - these might have duplicate
          triangles and so integration might be wrong

See also
    Foam::fieldValues
    Foam::functionObject

SourceFiles
    surfaceFieldValueBCFDK.C
    surfaceFieldValueBCFDKTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_surfaceFieldValueBCFDK_H
#define functionObjects_surfaceFieldValueBCFDK_H

#include "surfaceFieldValue.H"
#include "commonDictBCFDK.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace fieldValues
{

/*---------------------------------------------------------------------------*\
                      Class surfaceFieldValueBCFDK Declaration
\*---------------------------------------------------------------------------*/

class surfaceFieldValueBCFDK
:
    public surfaceFieldValue,
    public commonDictBCFDK
{

public:

    //- Run-time type information
    TypeName("surfaceFieldValueBCFDK");


    // Constructors

        //- Construct from name, Time and dictionary
        surfaceFieldValueBCFDK
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct from name, objectRegistry and dictionary
        surfaceFieldValueBCFDK
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );


    //- Destructor
    virtual ~surfaceFieldValueBCFDK();


    // Public Member Functions

        //- Custom templated helper function to output field values
        template<class Type>
        bool writeValuesBCFDK
        (
            const word& fieldName,
            const scalarField& weightField,
            const bool orient,
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
    #include "surfaceFieldValueBCFDKTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
