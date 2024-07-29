/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's createPatch.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2021-2021 FSD blueCAPE Lda  http://www.bluecape.com.pt
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
    reorderPatch

Description
    Utility to reorder patches of selected boundaries and their respective
    faces.

    It essentially uses the filterPatches function from createPatch.
\*---------------------------------------------------------------------------*/

#include "syncTools.H"
#include "argList.H"
#include "polyMesh.H"
#include "meshTools.H"
#include "IOdictionary.H"
#include "polyTopoChange.H"

#if defined( WIN32 ) || defined( WIN64 )
#include "Time.T.H"
#include "SortableList.T.H"
#include "IOPtrList.T.H"
#else
#include "Time.H"
#include "SortableList.H"
#include "IOPtrList.H"
#endif


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(IOPtrList<dictionary>, 0);
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

    const word oldInstance = mesh.pointsInstance();

    const word dictName("reorderPatchDict");
    #include "setSystemMeshDictionaryIO.H"

    Info<< "Reading " << dictName << nl << endl;

    IOdictionary dict(dictIO);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // If running parallel check same patches everywhere
    patches.checkParallelSync(true);

    polyTopoChange meshMod(mesh);

    // Read patch pairs from dictionary, for reordering
    List<Pair<word>> reorderPatchPairs
    (
        dict.lookup("reorderPatchPairs")
    );

    // 1. Compile list of patches being reordered
    wordList pNames(patches.names());
    labelList oldToNew(pNames.size());

    forAll(oldToNew, patchi)
    {
        oldToNew[patchi] = patchi;
    }

    Info<< "Lookup patch swaping indexes..." << nl<< endl;
    forAll(reorderPatchPairs, pairI)
    {
        label origin = findIndex(pNames, reorderPatchPairs[pairI].first());
        label target = findIndex(pNames, reorderPatchPairs[pairI].second());

        Info << "  Swapping '" << reorderPatchPairs[pairI].first() << "' for '"
             << reorderPatchPairs[pairI].second() << "', namely IDs " << origin
             << " with " << target << endl;

        label keep = oldToNew[target];
        oldToNew[target] = oldToNew[origin];
        oldToNew[origin] = keep;
    }

    // 2. Reorder patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Info<< "Reordering patches..." << nl<< endl;
    mesh.reorderPatches(oldToNew, false);

    Info<< "Doing topology modification to order faces." << nl << endl;
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, true);

    mesh.setInstance(oldInstance);

    // Write resulting mesh
    Info<< "Writing reordered boundary mesh to " << runTime.timeName() << nl << endl;
    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
