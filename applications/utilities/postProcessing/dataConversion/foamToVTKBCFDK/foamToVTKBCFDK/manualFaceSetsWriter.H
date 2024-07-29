/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
 2014-09-19 FSD blueCAPE Lda: Hack for writing faceSets along with allPatches
------------------------------------------------------------------------------
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
    Foam::manualFaceSetsWriter

Description
    Write patch fields

SourceFiles
    manualFaceSetsWriter.C

\*---------------------------------------------------------------------------*/

#ifndef manualFaceSetsWriter_H
#define manualFaceSetsWriter_H

#include "pointMesh.H"
#include "OFstream.H"
#include "volFields.H"
#include "pointFields.H"
#include "vtkMesh.H"
#include "indirectPrimitivePatch.H"

#if defined( WIN32 ) || defined( WIN64 )
#include "PrimitivePatchInterpolation.T.H"
#else
#include "PrimitivePatchInterpolation.H"
#endif

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class volPointInterpolation;

/*---------------------------------------------------------------------------*\
                           Class manualFaceSetsWriter Declaration
\*---------------------------------------------------------------------------*/

class manualFaceSetsWriter
{
    const vtkMesh& vMesh_;

    const bool binary_;

    const fileName fName_;

    labelList faceSetNameIDs_;
    labelList faceSetCellZoneIDs_;
    PtrList<primitiveFacePatch> faceSetsAsPrimitives_;

    std::ofstream os_;

    label nPoints_;

    label nFaces_;

public:

    // Constructors

        //- Construct from components
        manualFaceSetsWriter
        (
            const vtkMesh&,
            const bool binary,
            const fileName&,
            const PtrList<dictionary>& manuallyTaggedFaceSetsList
        );


    // Member Functions

        std::ofstream& os()
        {
            return os_;
        }

        label nPoints() const
        {
            return nPoints_;
        }

        label nFaces() const
        {
            return nFaces_;
        }

        //- Write nameIDs in nameID field
        void writeNameIDs();

        //- Write cellZoneIDs in cellZoneID field
        void writeCellZoneIDs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
