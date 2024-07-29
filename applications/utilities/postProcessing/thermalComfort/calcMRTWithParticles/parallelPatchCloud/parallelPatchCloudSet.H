/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's patchCloudSet.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2016-2018 FSD blueCAPE Lda  http://www.bluecape.com.pt
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
    Foam::parallelPatchCloudSet

Description
    Like cloudSet but samples nearest patch face, modified to account for a
    pure parallel run, namely where the positions are not global. This means
    that parallellization data has been stripped from the code, since it's not
    used.

SourceFiles
    parallelPatchCloudSet.C

\*---------------------------------------------------------------------------*/

#ifndef parallelPatchCloudSet_H
#define parallelPatchCloudSet_H

#include "sampledSet.H"

#if defined(WIN32) || defined(WIN64)
#include "DynamicList.T.H"
#include "HashSet.T.H"
#else
#include "DynamicList.H"
#include "HashSet.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class parallelPatchCloudSet Declaration
\*---------------------------------------------------------------------------*/

class parallelPatchCloudSet
:
    public sampledSet
{
    // Private Data

        //- Sampling points
        const List<point> sampleCoords_;

        //- Patches to sample
        const labelHashSet patchSet_;

        //- Maximum distance to look for nearest
        const scalar searchDist_;


    // Private Member Functions

        //- Samples all points in sampleCoords.
        void calcSamples
        (
            DynamicList<point>& samplingPts,
            DynamicList<label>& samplingCells,
            DynamicList<label>& samplingFaces,
            DynamicList<label>& samplingSegments,
            DynamicList<scalar>& samplingCurveDist
        ) const;

        //- Uses calcSamples to obtain samples. Copies them into *this.
        void genSamples();


public:

    //- Runtime type information
    TypeName("parallelPatchCloud");


    // Constructors

        //- Construct from components
        parallelPatchCloudSet
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const word& axis,
            const List<point>& sampleCoords,
            const labelHashSet& patchSet,
            const scalar searchDist
        );

        //- Construct from dictionary
        parallelPatchCloudSet
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const dictionary& dict
        );


    //- Destructor
    virtual ~parallelPatchCloudSet();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
