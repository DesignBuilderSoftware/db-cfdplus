/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Derived from OpenFOAM's snappyHexMesh library.
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
    Foam::snappySnapDriver

Description
    All to do with snapping to surface

SourceFiles
    snappySnapDriver.C
    snappySnapDriverFeature.C

\*---------------------------------------------------------------------------*/

#ifndef snappySnapDriver_H
#define snappySnapDriver_H

#include "meshRefinement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class motionSmoother;
class snapParameters;
class pointConstraint;

/*---------------------------------------------------------------------------*\
                           Class snappySnapDriver Declaration
\*---------------------------------------------------------------------------*/

class snappySnapDriver
{
    // Private Data

        //- Mesh+surface
        meshRefinement& meshRefiner_;

        //- From global surface region to master side patch
        const labelList globalToMasterPatch_;

        //- From global surface region to slave side patch
        const labelList globalToSlavePatch_;


    // Private Member Functions


        // Snapping

            //- Calculates (geometric) shared points
            //  Requires PackedList to be sized and initialised
            static label getCollocatedPoints
            (
                const scalar tol,
                const pointField&,
                PackedBoolList&
            );

            //- Calculate displacement per patch point to smooth out patch.
            //  Quite complicated in determining which points to move where.
            static pointField smoothPatchDisplacement
            (
                const motionSmoother&,
                const List<labelPair>&
            );

            //- Check that face zones are synced
            void checkCoupledFaceZones() const;

            //- Per edge distance to patch
            static tmp<scalarField> edgePatchDist
            (
                const pointMesh&,
                const indirectPrimitivePatch&
            );

            //- Write displacement as .obj file.
            static void dumpMove
            (
                const fileName&,
                const pointField&,
                const pointField&
            );

            //- Check displacement is outwards pointing
            static bool outwardsDisplacement
            (
                const indirectPrimitivePatch&,
                const vectorField&
            );

            //- Detect warpage
            void detectWarpedFaces
            (
                const scalar featureCos,
                const indirectPrimitivePatch& pp,

                DynamicList<label>& splitFaces,
                DynamicList<labelPair>& splits
            ) const;

            // Feature line snapping

                //- Is point on two feature edges that make a largish angle?
                bool isFeaturePoint
                (
                    const scalar featureCos,
                    const indirectPrimitivePatch& pp,
                    const PackedBoolList& isFeatureEdge,
                    const label pointi
                ) const;

                void smoothAndConstrain
                (
                    const PackedBoolList& isMasterEdge,
                    const indirectPrimitivePatch& pp,
                    const labelList& meshEdges,
                    const List<pointConstraint>& constraints,
                    vectorField& disp
                ) const;

                void calcNearest
                (
                    const label iter,
                    const indirectPrimitivePatch& pp,
                    vectorField& pointDisp,
                    vectorField& pointSurfaceNormal,
                    vectorField& pointRotation
                ) const;

                void calcNearestFace
                (
                    const label iter,
                    const indirectPrimitivePatch& pp,
                    const scalarField& faceSnapDist,
                    vectorField& faceDisp,
                    vectorField& faceSurfaceNormal,
                    labelList& faceSurfaceRegion,
                    vectorField& faceRotation
                ) const;

                //- Collect (possibly remote) per point data of all surrounding
                //  faces
                //  - faceSurfaceNormal
                //  - faceDisp
                //  - faceCentres&faceNormal
                void calcNearestFacePointProperties
                (
                    const label iter,
                    const indirectPrimitivePatch& pp,

                    const vectorField& faceDisp,
                    const vectorField& faceSurfaceNormal,
                    const labelList& faceSurfaceRegion,

                    List<List<point>>& pointFaceSurfNormals,
                    List<List<point>>& pointFaceDisp,
                    List<List<point>>& pointFaceCentres,
                    List<labelList>&    pointFacePatchID
                ) const;

                //- Gets passed in offset to nearest point on feature
                //  edge. Calculates if the point has a different number of
                //  faces on either side of the feature and if so attracts the
                //  point to that non-dominant plane.
                void correctAttraction
                (
                    const DynamicList<point>& surfacePoints,
                    const DynamicList<label>& surfaceCounts,
                    const point& edgePt,
                    const vector& edgeNormal,   // normalised normal
                    const point& pt,
                    vector& edgeOffset  // offset from pt to point on edge
                ) const;


                //- For any reverse (so from feature back to mesh) attraction:
                //  add attraction if diagonal points on face attracted
                void stringFeatureEdges
                (
                    const label iter,
                    const scalar featureCos,

                    const indirectPrimitivePatch& pp,
                    const scalarField& snapDist,

                    const vectorField& rawPatchAttraction,
                    const List<pointConstraint>& rawPatchConstraints,

                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                //- Remove constraints of points next to multi-patch points
                //  to give a bit more freedom of the mesh to conform to the
                //  multi-patch points. Bit dodgy for simple cases.
                void releasePointsNextToMultiPatch
                (
                    const label iter,
                    const scalar featureCos,

                    const indirectPrimitivePatch& pp,
                    const scalarField& snapDist,

                    const List<List<point>>& pointFaceCentres,
                    const labelListList& pointFacePatchID,

                    const vectorField& rawPatchAttraction,
                    const List<pointConstraint>& rawPatchConstraints,

                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                //- Detect any diagonal attraction. Returns indices in face
                //  or (-1, -1) if none
                labelPair findDiagonalAttraction
                (
                    const indirectPrimitivePatch& pp,
                    const vectorField& patchAttraction,
                    const List<pointConstraint>& patchConstraints,
                    const label facei
                ) const;

                //- Avoid attraction across face diagonal since would
                //  cause face squeeze
                void avoidDiagonalAttraction
                (
                    const label iter,
                    const scalar featureCos,
                    const indirectPrimitivePatch& pp,
                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                //- Return hit if on multiple points
                pointIndexHit findMultiPatchPoint
                (
                    const point& pt,
                    const labelList& patchIDs,
                    const List<point>& faceCentres
                ) const;

                //- Return hit if faces-on-the-same-normalplane are on multiple
                //  patches
                //  - false, index=-1 : single patch
                //  - true , index=0  : multiple patches but on different
                //                      normals planes (so geometric feature
                //                      edge is also a region edge)
                //  - true , index=1  : multiple patches on same normals plane
                //                      i.e. flat region edge
                pointIndexHit findMultiPatchPoint
                (
                    const point& pt,
                    const labelList& pfPatchID,
                    const DynamicList<vector>& surfaceNormals,
                    const labelList& faceToNormalBin
                ) const;

                //- Return index of similar normal
                label findNormal
                (
                    const scalar featureCos,
                    const vector& faceSurfaceNormal,
                    const DynamicList<vector>& surfaceNormals
                ) const;

                //- Determine attraction and constraints for single point
                //  using sampled surrounding of the point
                void featureAttractionUsingReconstruction
                (
                    const label iter,
                    const scalar featureCos,

                    const indirectPrimitivePatch& pp,
                    const scalarField& snapDist,
                    const vectorField& nearestDisp,
                    const label pointi,

                    const List<List<point>>& pointFaceSurfNormals,
                    const List<List<point>>& pointFaceDisp,
                    const List<List<point>>& pointFaceCentres,
                    const labelListList& pointFacePatchID,

                    DynamicList<point>& surfacePoints,
                    DynamicList<vector>& surfaceNormals,
                    labelList& faceToNormalBin,

                    vector& patchAttraction,
                    pointConstraint& patchConstraint
                ) const;

                //- Determine attraction and constraints for all points
                //  using sampled surrounding of the point
                void featureAttractionUsingReconstruction
                (
                    const label iter,
                    const bool avoidSnapProblems,
                    const scalar featureCos,
                    const indirectPrimitivePatch& pp,
                    const scalarField& snapDist,
                    const vectorField& nearestDisp,

                    const List<List<point>>& pointFaceSurfNormals,
                    const List<List<point>>& pointFaceDisp,
                    const List<List<point>>& pointFaceCentres,
                    const labelListList& pointFacePatchID,

                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                //- Determine geometric features and attraction to equivalent
                //  surface features
                void determineFeatures
                (
                    const label iter,
                    const scalar featureCos,
                    const bool multiRegionFeatureSnap,

                    const indirectPrimitivePatch&,
                    const scalarField& snapDist,
                    const vectorField& nearestDisp,

                    const List<List<point>>& pointFaceSurfNormals,
                    const List<List<point>>& pointFaceDisp,
                    const List<List<point>>& pointFaceCentres,
                    const labelListList& pointFacePatchID,

                    List<labelList>& pointAttractor,
                    List<List<pointConstraint>>& pointConstraints,
                    // Feature-edge to pp point
                    List<List<DynamicList<point>>>& edgeAttractors,
                    List<List<DynamicList<pointConstraint>>>& edgeConstraints,
                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                //- Determine features originating from bafles and
                //  and add attraction to equivalent surface features
                void determineBaffleFeatures
                (
                    const label iter,
                    const scalar featureCos,

                    const indirectPrimitivePatch& pp,
                    const scalarField& snapDist,

                    // Feature-point to pp point
                    List<labelList>& pointAttractor,
                    List<List<pointConstraint>>& pointConstraints,
                    // Feature-edge to pp point
                    List<List<DynamicList<point>>>& edgeAttractors,
                    List<List<DynamicList<pointConstraint>>>& edgeConstraints,
                    // pp point to nearest feature
                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                //- Patch corner readjustments for those in the same plane
                //  surface. Hack put in by FSD blueCAPE Lda on June 2014,
                //  later revised in April 2018.
                void patchCornerAdjustments
                (
                    const bool multiRegionFeatureSnap,

                    const indirectPrimitivePatch& pp,
                    const scalarField& snapDist,

                    const refinementFeatures& features,

                    // Feature-point to pp point
                    List<labelList>& pointAttractor,
                    List<List<pointConstraint>>& pointConstraints
                ) const;

                void reverseAttractMeshPoints
                (
                    const label iter,

                    const indirectPrimitivePatch& pp,
                    const scalarField& snapDist,

                    // Feature-point to pp point
                    const List<labelList>& pointAttractor,
                    const List<List<pointConstraint>>& pointConstraints,
                    // Feature-edge to pp point
                    const List<List<DynamicList<point>>>& edgeAttractors,
                    const List<List<DynamicList<pointConstraint>>>&,

                    const vectorField& rawPatchAttraction,
                    const List<pointConstraint>& rawPatchConstraints,

                    // pp point to nearest feature
                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                //- Find point on nearest feature edge (within searchDist).
                //  Return point and feature
                //  and store feature-edge to mesh-point and vice versa
                Tuple2<label, pointIndexHit> findNearFeatureEdge
                (
                    const bool isRegionEdge,

                    const indirectPrimitivePatch& pp,
                    const scalarField& snapDist,
                    const label pointi,
                    const point& estimatedPt,

                    List<List<DynamicList<point>>>&,
                    List<List<DynamicList<pointConstraint>>>&,
                    vectorField&,
                    List<pointConstraint>&
                ) const;

                //- Find nearest feature point (within searchDist).
                //  Return feature point
                //  and store feature-point to mesh-point and vice versa.
                //  If another mesh point already referring to this feature
                //  point and further away, reset that one to a near feature
                //  edge (using findNearFeatureEdge above)
                Tuple2<label, pointIndexHit> findNearFeaturePoint
                (
                    const bool isRegionEdge,

                    const indirectPrimitivePatch& pp,
                    const scalarField& snapDist,
                    const label pointi,
                    const point& estimatedPt,

                    // Feature-point to pp point
                    List<labelList>& pointAttractor,
                    List<List<pointConstraint>>& pointConstraints,
                    // Feature-edge to pp point
                    List<List<DynamicList<point>>>& edgeAttractors,
                    List<List<DynamicList<pointConstraint>>>& edgeConstraints,
                    // pp point to nearest feature
                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                void featureAttractionUsingFeatureEdges
                (
                    const label iter,
                    const bool avoidSnapProblems,
                    const scalar featureCos,
                    const bool multiRegionFeatureSnap,
                    const indirectPrimitivePatch& pp,
                    const scalarField& snapDist,
                    const vectorField& nearestDisp,

                    const List<List<point>>& pointFaceSurfNormals,
                    const List<List<point>>& pointFaceDisp,
                    const List<List<point>>& pointFaceCentres,
                    const labelListList& pointFacePatchID,

                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                void preventFaceSqueeze
                (
                    const label iter,
                    const scalar featureCos,
                    const indirectPrimitivePatch& pp,
                    const scalarField& snapDist,

                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                //- Top level feature attraction routine. Gets given
                //  displacement to nearest surface in nearestDisp
                //  and calculates new displacement taking into account
                //  features
                vectorField calcNearestSurfaceFeature
                (
                    const snapParameters& snapParams,
                    const bool avoidSnapProblems,
                    const label iter,
                    const scalar featureCos,
                    const scalar featureAttract,
                    const scalarField& snapDist,
                    const vectorField& nearestDisp,
                    motionSmoother& meshMover,
                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                static void writeUnsnappedToVTK
                (
                    const fileName& vtkFilePath,
                    const label vtkIter,
                    const DynamicList<vector>& pointList,
                    const DynamicList<scalar>& minNonSnappedDist
                );

public:

    //- Runtime type information
    ClassName("snappySnapDriver");


    // Constructors

        //- Construct from components
        snappySnapDriver
        (
            meshRefinement& meshRefiner,
            const labelList& globalToMasterPatch,
            const labelList& globalToSlavePatch
        );

        //- Disallow default bitwise copy construction
        snappySnapDriver(const snappySnapDriver&) = delete;


    // Member Functions

        // Snapping

            //- Merge baffles.
            autoPtr<mapPolyMesh> mergeZoneBaffles(const List<labelPair>&);

            //- Calculate edge length per patch point.
            static scalarField calcSnapDistance
            (
                const fvMesh& mesh,
                const snapParameters& snapParams,
                const indirectPrimitivePatch&
            );

            //- Smooth the mesh (patch and internal) to increase visibility
            //  of surface points (on castellated mesh) w.r.t. surface.
            static void preSmoothPatch
            (
                const meshRefinement& meshRefiner,
                const snapParameters& snapParams,
                const label nInitErrors,
                const List<labelPair>& baffles,
                motionSmoother&
            );

            //- Get points both on patch and facezone.
            static labelList getZoneSurfacePoints
            (
                const fvMesh& mesh,
                const indirectPrimitivePatch&,
                const word& zoneName
            );

            //- Helper: calculate average cell centre per point
            static tmp<pointField> avgCellCentres
            (
                const fvMesh& mesh,
                const indirectPrimitivePatch&
            );

            //- Per patch point override displacement if in gap situation
            void detectNearSurfaces
            (
                const scalar planarCos,
                const indirectPrimitivePatch&,
                const pointField& nearestPoint,
                const vectorField& nearestNormal,
                vectorField& disp
            ) const;

            //- Per patch point calculate point on nearest surface. Set as
            //  boundary conditions of motionSmoother displacement field. Return
            //  displacement of patch points.
            static vectorField calcNearestSurface
            (
                const meshRefinement& meshRefiner,
                const scalarField& snapDist,
                const indirectPrimitivePatch&,
                pointField& nearestPoint,
                vectorField& nearestNormal,
                const fileName &vtkFilePath = fileName(""),
                const label vtkIter = 0
            );

            //- Smooth the displacement field to the internal.
            void smoothDisplacement
            (
                const snapParameters& snapParams,
                motionSmoother&
            ) const;

            //- Do the hard work: move the mesh according to displacement,
            //  locally relax the displacement. Return true if ended up with
            //  correct mesh, false if not.
            bool scaleMesh
            (
                const snapParameters& snapParams,
                const label nInitErrors,
                const List<labelPair>& baffles,
                motionSmoother&
            );

            //- Repatch faces according to surface nearest the face centre
            //  - calculate face-wise snap distance as max of point-wise
            //  - calculate face-wise nearest surface point
            //  - repatch face according to patch for surface point.
            autoPtr<mapPolyMesh> repatchToSurface
            (
                const snapParameters& snapParams,
                const labelList& adaptPatchIDs,
                const labelList& preserveFaces
            );

            void doSnap
            (
                const dictionary& snapDict,
                const dictionary& motionDict,
                const scalar featureCos,
                const scalar planarAngle,
                const snapParameters& snapParams
            );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const snappySnapDriver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
