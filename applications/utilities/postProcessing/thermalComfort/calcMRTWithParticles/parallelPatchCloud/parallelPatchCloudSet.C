/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "parallelPatchCloudSet.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "treeBoundBox.H"
#include "treeDataFace.H"
#if defined(WIN32) || defined(WIN64)
#include "Time.T.H"
#else
#include "Time.H"
#endif
#include "meshTools.H"
// For 'nearInfo' helper class only
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(parallelPatchCloudSet, 0);
    addToRunTimeSelectionTable(sampledSet, parallelPatchCloudSet, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::parallelPatchCloudSet::calcSamples
(
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
) const
{
    if (debug)
    {
        Pout<< "parallelPatchCloudSet : sampling on patches :" << endl;
    }

    // Construct search tree for all patch faces.
    label sz = 0;
    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        const polyPatch& pp = mesh().boundaryMesh()[iter.key()];

        sz += pp.size();

        if (debug)
        {
            Pout<< "    " << pp.name() << " size " << pp.size() << endl;
        }
    }

    labelList patchFaces(sz);
    sz = 0;
    treeBoundBox bb(point::max, point::min);
    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        const polyPatch& pp = mesh().boundaryMesh()[iter.key()];

        forAll(pp, i)
        {
            patchFaces[sz++] = pp.start()+i;
        }

        // Do not do reduction.
        const boundBox patchBb(pp.points(), pp.meshPoints(), false);

        bb.min() = min(bb.min(), patchBb.min());
        bb.max() = max(bb.max(), patchBb.max());
    }

    // Make bb asymetric just to avoid problems on symmetric meshes
    bb = bb.extend(1e-4);

    // Make sure bb is 3D.
    bb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
    bb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);


    indexedOctree<treeDataFace> patchTree
    (
        treeDataFace    // all information needed to search faces
        (
            false,      // do not cache bb
            mesh(),
            patchFaces  // boundary faces only
        ),
        bb,             // overall search domain
        8,              // maxLevel
        10,             // leafsize
        3.0             // duplicity
    );

    // Force calculation of face-diagonal decomposition
    (void)mesh().tetBasePtIs();


    // All the info for nearest. Construct to miss
    List<pointIndexHit> nearest(sampleCoords_.size());

    forAll(sampleCoords_, sampleI)
    {
        const point& sample = sampleCoords_[sampleI];

        pointIndexHit& nearInfo = nearest[sampleI];

        // Find the nearest locally
        if (patchFaces.size())
        {
            nearInfo = patchTree.findNearest(sample, sqr(searchDist_));
        }
        else
        {
            nearInfo.setMiss();
        }


        // Fill in the distance field and the processor field
        if (nearInfo.hit())
        {
            // Set nearest to mesh face label
            nearInfo.setIndex(patchFaces[nearInfo.index()]);
        }
    }


    if (debug)
    {
        OFstream str
        (
            mesh().time().path()
          / name()
          + "_nearest.obj"
        );
        Info<< "Dumping mapping as lines from supplied points to"
            << " nearest patch face to file " << str.name() << endl;

        label vertI = 0;

        forAll(nearest, i)
        {
            if (nearest[i].hit())
            {
                meshTools::writeOBJ(str, sampleCoords_[i]);
                vertI++;
                meshTools::writeOBJ(str, nearest[i].hitPoint());
                vertI++;
                str << "l " << vertI-1 << ' ' << vertI << nl;
            }
        }
    }


    // Store the sampling locations on the nearest processor
    forAll(nearest, sampleI)
    {
        const pointIndexHit& nearInfo = nearest[sampleI];

        if (nearInfo.hit())
        {
            label facei = nearInfo.index();

            samplingPts.append(nearInfo.hitPoint());
            samplingCells.append(mesh().faceOwner()[facei]);
            samplingFaces.append(facei);
            samplingSegments.append(0);
            samplingCurveDist.append(1.0 * sampleI);
        }
        else
        {
            // Position not found. Mark with special value
            // which is intercepted when interpolating
            samplingPts.append(sampleCoords_[sampleI]);
            samplingCells.append(-1);
            samplingFaces.append(-1);
            samplingSegments.append(0);
            samplingCurveDist.append(1.0 * sampleI);
        }
    }
}


void Foam::parallelPatchCloudSet::genSamples()
{
    // Storage for sample points
    DynamicList<point> samplingPts;
    DynamicList<label> samplingCells;
    DynamicList<label> samplingFaces;
    DynamicList<label> samplingSegments;
    DynamicList<scalar> samplingCurveDist;

    calcSamples
    (
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );

    samplingPts.shrink();
    samplingCells.shrink();
    samplingFaces.shrink();
    samplingSegments.shrink();
    samplingCurveDist.shrink();

    setSamples
    (
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parallelPatchCloudSet::parallelPatchCloudSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const word& axis,
    const List<point>& sampleCoords,
    const labelHashSet& patchSet,
    const scalar searchDist
)
:
    sampledSet(name, mesh, searchEngine, axis),
    sampleCoords_(sampleCoords),
    patchSet_(patchSet),
    searchDist_(searchDist)
{
    genSamples();

    if (debug)
    {
        write(Pout);
    }
}


Foam::parallelPatchCloudSet::parallelPatchCloudSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    sampleCoords_(dict.lookup("points")),
    patchSet_
    (
        mesh.boundaryMesh().patchSet
        (
            wordReList(dict.lookup("patches"))
        )
    ),
    searchDist_(dict.lookup<scalar>("maxDistance"))
{
    genSamples();

    if (debug)
    {
        write(Pout);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::parallelPatchCloudSet::~parallelPatchCloudSet()
{}


// ************************************************************************* //
