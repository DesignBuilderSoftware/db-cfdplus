/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on blueCFD(R)-Kernel's reorderPatch.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2023-2023 FSD blueCAPE Lda  http://www.bluecape.com.pt
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

Application
    checkPatchStats

Description
    Utility to check selected stats for patches and faceSets, regarding areas,
    centroids and normals.

    Used reorderPatch and createBafflesBCFDK as a basis.
\*---------------------------------------------------------------------------*/

#include "syncTools.H"
#include "argList.H"
#include "polyMesh.H"
#include "meshTools.H"
#include "IOdictionary.H"
#include "faceSet.H"

#if defined( WIN32 ) || defined( WIN64 )
#include "Time.T.H"
#else
#include "Time.H"
#endif


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

enum surfaceType
{
    patchEnum,
    faceSetEnum
};

const NamedEnum<surfaceType, 2> surfaceTypeNames;

template<> const char* Foam::NamedEnum<surfaceType, 2>::names[] =
{
    "patch",
    "faceSet"
};

template<class patchObj> void processPatchObject
(
    const patchObj& patch,
    const word& surfaceName,
    const vector& surfCentre,
    const vector& surfNormal
)
{
    const vectorField faceCentres(patch.faceCentres());
    const vectorField areasSf(patch.faceAreas());
    const scalarField areas(mag(patch.faceAreas()));
    const scalar totalArea = sum(areas);

    vector centroid = sum(faceCentres*areas);
    if (totalArea > vSmall)
    {
        centroid /= totalArea;
    }

    vectorField normals(areasSf/areas);

    // Given that geometry surface normals don't always respect the
    // patch/faceSet normals, we have to take into account the one with the
    // smallest error; and also have to account for the fact that faceSets
    // don't always respect the face orientation of the surface that gave it
    // shape, so we need to apply a similar correction as done in
    // orientFaceZones.

    const vector outsideSurfPoint = surfNormal + surfCentre;

    forAll (normals, facei)
    {

        const vector d = outsideSurfPoint - faceCentres[facei];
        vector& fn = normals[facei];

        if ((fn&d) < 0)
        {
            normals[facei] *= scalar(-1.0);
        }
    }

    vectorField diffNormals(normals - surfNormal);
    scalar cumulativeDifference = sumMag(diffNormals*areas);

    vector centreDiff = centroid - surfCentre;

    Info
        << surfaceName <<  token::COLON
        << token::SPACE
        << "area" << token::ASSIGN << totalArea << token::END_STATEMENT
        << token::SPACE
        << "sumAreaNDiff" << token::ASSIGN << cumulativeDifference
            << token::END_STATEMENT
        << token::SPACE
        << "centroidDiff" << token::ASSIGN << centreDiff
            << token::END_STATEMENT
        << endl;

    // Debug normals:
    //Info << nl
    //    << "normals: " << normals << endl;
}

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    #include "addDictOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();

    Foam::word meshRegionName = polyMesh::defaultRegion;
    args.optionReadIfPresent("region", meshRegionName);

    #include "createNamedPolyMesh.H"

    const word dictName("checkPatchStatsDict");
    #include "setSystemMeshDictionaryIO.H"

    Info<< "Reading " << dictName << nl << endl;

    IOdictionary dict(dictIO);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // If running parallel check same patches everywhere
    patches.checkParallelSync(true);

    Info << "Surface stats:" << endl;

    forAllConstIter(dictionary, dict, item)
    {
        if (item().isDict())
        {
            // Process this surface
            const word& surfaceName = item().keyword();
            const dictionary& itemDict = item().dict();

            const surfaceType surfType = surfaceTypeNames
            [
                itemDict.lookup<word>("type")
            ];

            const vector surfCentre(itemDict.lookup<vector>("centre"));
            const vector surfNormal(itemDict.lookup<vector>("normal"));

            if (surfType == patchEnum)
            {
                const polyPatch& patch = patches[surfaceName];

                if(patch.size()>0)
                {
                    processPatchObject
                    (
                        patch,
                        surfaceName,
                        surfCentre,
                        surfNormal
                    );
                }
            }
            else if (surfType == faceSetEnum)
            {
                const fileName pathToSet
                (
                    mesh.time().path()/topoSet::localPath(mesh, surfaceName)
                );

                if(exists(pathToSet))
                {
                    faceSet set(mesh, surfaceName);

                    if (set.size()>0)
                    {

                        const indirectPrimitivePatch setPatch
                        (
                            IndirectList<face>(mesh.faces(), set.sortedToc()),
                            mesh.points()
                        );

                        processPatchObject
                        (
                            setPatch,
                            surfaceName,
                            surfCentre,
                            surfNormal
                        );
                    }
                }
            }

        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
