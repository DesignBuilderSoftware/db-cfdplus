/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Derived from OpenFOAM's surfaceFieldValue and probes function object.
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
    Foam::functionObjects::interpolatingProbesBCFDK

Group
    grpFieldFunctionObjects

Description
    Provides a single probe position, including interpolation scheme.

    Examples of function object specification:
    \verbatim
    probe1
    {
        type            interpolatingProbesBCFDK;
        libs            ("libMonitorsBCFDK.so");

        writeControl    writeTime;

        interpolationScheme  cellPoint;

        probePosition  (0.0254 0.0253 0.0);

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
        type         | type name: interpolatingProbesBCFDK   | yes         |
        name         | name of face regionType if required  | no |
        probePosition | vector with the position of the probe | yes |
        fields       | list of fields to operate on | yes    |
        interpolationScheme| list of fields to operate on | yes    |
        identifiers  | corresponding list of names for each field entry | yes |
        modes        | list of calculation modes:
                       scalar, magnitude, componentX,
                       componentY, componentZ | yes | magnitude
    \endtable

    Choice of \c interpolationScheme:
      - cell          : use cell-centre value only; constant over cells
                           (default)
      - cellPoint     : use cell-centre and vertex values
      - cellPointFace : use cell-centre, vertex and face values.
      - pointMVC      : use point values only (Mean Value Coordinates)
      - cellPatchConstrained : like 'cell' but uses cell-centre except on
                               boundary faces where it uses the boundary value.
                               For use with e.g. patchCloudSet.
    1. vertex values determined from neighbouring cell-centre values
    2. face values determined using the current face interpolation scheme
       for the field (linear, gamma, etc.)

See also
    Foam::functionObject

SourceFiles
    interpolatingProbesBCFDK.C
    interpolatingProbesBCFDKTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_interpolatingProbesBCFDK_H
#define functionObjects_interpolatingProbesBCFDK_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "commonDictBCFDK.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                      Class interpolatingProbesBCFDK Declaration
\*---------------------------------------------------------------------------*/

class interpolatingProbesBCFDK
:
    public fvMeshFunctionObject,
    public logFiles,
    public commonDictBCFDK
{

protected:

    // Protected data

        //- List of fields
        wordList fields_;

        //- Interpolation scheme to use
        word interpolationScheme_;

        //- Position of the probe
        vector probePosition_;

        // Cell to be probed (obtained from the probe position)
        label elementProbe_;

        // Face to be probed
        label faceProbe_;

    // Protected Member Functions

        //- Sample and report the point position
        template<class Type>
        bool probePosition
        (
            const word& fieldName,
            const modeType& mode,
            const word& identifier
        );

        //- Find cells and faces containing probes
        //  Copied from probes::findElements, because it was too cumbersome to
        //  inherit the whole class.
        virtual void findElement(const fvMesh&);

        //- Write file header, implementing from logFiles
        void writeFileHeader(const label i);

public:

    //- Run-time type information
    TypeName("interpolatingProbesBCFDK");


    // Constructors

        //- Construct from name, Time and dictionary
        interpolatingProbesBCFDK
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct from name, objectRegistry and dictionary
        interpolatingProbesBCFDK
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );


    //- Destructor
    virtual ~interpolatingProbesBCFDK();


    // Public Member Functions

        //- Read from dictionary
        virtual bool read(const dictionary&);

        //- Do nothing on execute
        virtual bool execute();

        //- Calculate and write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "interpolatingProbeBCFDKTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
