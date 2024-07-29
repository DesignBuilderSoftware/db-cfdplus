
//FIXME: "searchDist" should be an option in the Dict
const scalar searchDist(0.0001);

meshSearch searchEngine(mesh);

List<point> sampleCoords(kinematicCloud.size());

label counter = 0;
forAllIter(basicKinematicCollidingCloud, kinematicCloud, pIter)
{
    basicKinematicCollidingCloud::parcelType& p = pIter();

    sampleCoords[counter] = p.position();
    counter++;
}

// Sanity check for total existing particles versus total number of cells
{
    counter = returnReduce(counter, sumOp<label>());

    label cellCount = mesh.nCells();
    cellCount = returnReduce(cellCount, sumOp<label>());

    if (cellCount != counter)
    {
      Info<< "WARNING: Several particles ignored due to flaws in the mesh."
        << " Total number: " << (cellCount - counter)
        << endl;
    }
}

parallelPatchCloudSet samplingCloud
(
    kinematicCloudName,
    mesh,
    searchEngine,
    "xyz",
    sampleCoords,
    mesh.boundaryMesh().patchSet
    (
        wordReList(1,wordRe(".*",wordRe::compOption::regExp)) //all patches
    ),
    searchDist
);

// Relevant options (from OpenFOAM's sampleDict):
//
//    cellPointFace : use cell-centre, vertex and face values.
//    pointMVC      : use point values only (Mean Value Coordinates)
//    cellPatchConstrained : like 'cell' but uses cell-centre except on
//                           boundary faces where it uses the boundary value.
//                           For use with e.g. patchCloudSet.
// Hard-coded, since we only want patch-based values.
word interpolationScheme("cellPatchConstrained");

autoPtr<interpolation<scalar> > interpolator
(
    interpolation<scalar>::New(interpolationScheme, T)
);

scalarField sampledT(samplingCloud.size());
forAll(samplingCloud, sampleI)
{
    const point& samplePt = samplingCloud[sampleI];
    label cellI = samplingCloud.cells()[sampleI];
    label faceI = samplingCloud.faces()[sampleI];

    if (cellI == -1 && faceI == -1)
    {
        // Special condition for illegal sampling points
        sampledT[sampleI] = 0.0;
    }
    else
    {
        sampledT[sampleI] = interpolator().interpolate
        (
            samplePt,
            cellI,
            faceI
        );
    }
}
