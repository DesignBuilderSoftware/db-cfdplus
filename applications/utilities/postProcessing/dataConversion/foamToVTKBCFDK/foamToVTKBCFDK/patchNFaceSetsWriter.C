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

Modifications
    This file has been created by blueCAPE's unofficial mingw patches for
    OpenFOAM.
    For more information about these patches, visit:
        http://www.bluecape.com.pt/blueCFD

    Modifications made:
      - Derived from the patches for blueCFD 2.1, which in turn were derived
        from 2.0 and 1.7.
      - Always open the files in binary mode, because of how things work on
        Windows.
        - Note: This modification is hard to stipulate who implemented this
                first. Symscape's port only began integrating this fix after
                blueCFD 1.7-2 was released with this sort of specific fix.
                But it was Symscape that first implemented this kind of
                "binary" fix in OpenFOAM's core code.

\*---------------------------------------------------------------------------*/

#include "patchNFaceSetsWriter.H"
#include "vtkWriteFieldOps.H"
#include "faceSet.H"
#include "stringListOps.H"

#if defined( WIN32 ) || defined( WIN64 )
#include "IOobjectList.T.H"
#else
#include "IOobjectList.H"
#endif

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchNFaceSetsWriter::patchNFaceSetsWriter
(
    const vtkMesh& vMesh,
    const bool binary,
    const bool nearCellValue,
    const fileName& fName,
    const labelList& patchIDs,
    const List<wordRe>& faceSetsToExclude,
    ListHashTable<vector,word> faceSetNormals
)
:
    vMesh_(vMesh),
    binary_(binary),
    nearCellValue_(nearCellValue),
    fName_(fName),
    patchIDs_(patchIDs),
    faceSetIDs_(),
    faceSetsAsPrimitives_(),
    os_(fName.c_str(),
        std::ios_base::out|std::ios_base::binary) //a must for Windows!
{
    const fvMesh& mesh = vMesh_.mesh();
    const faceList& meshFaces = mesh.faces();
    const pointField& meshPoints = mesh.points();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    PtrList<const faceSet> faceSets;

    Info<< "Exporting combined patches and faceSets (renumbered IDs):" << endl;

    forAll(patchIDs_, i)
    {
        const polyPatch& pp = patches[patchIDs_[i]];

        Info<< "    patch " << i << " " << pp.name() << endl;
    }

    // Read sets
    // Note: All sets are read, but omitted sets are tagged with ID = -1
    IOobjectList objects(mesh.thisDb(), mesh.facesInstance(), "polyMesh/sets");
    {
        label lastID = patchIDs_.size();
        label faceSetI = 0;

        IOobjectList fSets(objects.lookupClass(faceSet::typeName));

        faceSetIDs_.setSize(fSets.size());
        faceSets.setSize(fSets.size());
        faceSetsAsPrimitives_.setSize(fSets.size());

        forAllConstIter(IOobjectList, fSets, iter)
        {
            faceSets.set(faceSetI, new faceSet(*iter()));

            if(!findStrings(faceSetsToExclude, faceSets[faceSetI].name()))
            {
                faceSetIDs_[faceSetI] = lastID;

                Info<< "    faceSet " << lastID
                    << " " << faceSets[faceSetI].name() << endl;
                lastID++;
            }
            else
            {
                faceSetIDs_[faceSetI] = -1;
            }

            faceSetI++;
        }
    }

    // Write header
    vtkWriteOps::writeHeader(os_, binary_, "patches");
    os_ << "DATASET POLYDATA" << std::endl;

    // Write topology
    nPoints_ = 0;
    nFaces_ = 0;
    label nFaceVerts = 0;

    forAll(patchIDs_, i)
    {
        const polyPatch& pp = patches[patchIDs_[i]];

        nPoints_ += pp.nPoints();
        nFaces_ += pp.size();

        forAll(pp, faceI)
        {
            nFaceVerts += pp[faceI].size() + 1;
        }
    }

    forAll(faceSets, set)
    {
        const faceSet &fs = faceSets[set];

        if(faceSetIDs_[set] != -1)
        {
            faceList setFaces(fs.size());
            label setFaceI = 0;

            forAllConstIter(faceSet, fs, iter)
            {
                setFaces[setFaceI] = meshFaces[iter.key()];
                setFaceI++;
            }

            faceSetsAsPrimitives_.set(
                set,
                new primitiveFacePatch(setFaces, meshPoints)
                );

            const primitiveFacePatch& pp = faceSetsAsPrimitives_[set];

            nPoints_ += pp.nPoints();
            nFaces_ += pp.size();

            forAll(pp, faceI)
            {
                nFaceVerts += pp[faceI].size() + 1;
            }
        }
    }

    os_ << "POINTS " << nPoints_ << " float" << std::endl;

    DynamicList<floatScalar> ptField(3*nPoints_);

    forAll(patchIDs_, i)
    {
        const polyPatch& pp = patches[patchIDs_[i]];

        vtkWriteOps::insert(pp.localPoints(), ptField);
    }

    forAll(faceSetsAsPrimitives_, i)
    {
        if(faceSetIDs_[i] != -1)
        {
            const primitiveFacePatch& pp = faceSetsAsPrimitives_[i];

            vtkWriteOps::insert(pp.localPoints(), ptField);
        }
    }

    vtkWriteOps::write(os_, binary_, ptField);

    os_ << "POLYGONS " << nFaces_ << ' ' << nFaceVerts << std::endl;

    DynamicList<label> vertLabels(nFaceVerts);

    label offset = 0;

    forAll(patchIDs_, i)
    {
        const polyPatch& pp = patches[patchIDs_[i]];

        forAll(pp, faceI)
        {
            const face& f = pp.localFaces()[faceI];

            vertLabels.append(f.size());
            vtkWriteOps::insert(f + offset, vertLabels);
        }
        offset += pp.nPoints();
    }

    forAll(faceSetsAsPrimitives_, i)
    {
        if(faceSetIDs_[i] != -1)
        {
            const primitiveFacePatch& pp = faceSetsAsPrimitives_[i];

            ListHashTable<vector,word>::const_iterator iteFSN =
                faceSetNormals.find(faceSets[i].name());
            const bool normalFound = iteFSN != faceSetNormals.end();

            if(normalFound)
            {
                const vector& surfaceNormal(*iteFSN);

                forAll(pp, faceI)
                {
                    const face& f = pp.localFaces()[faceI];
                    const vector& faceNormal = pp.faceNormals()[faceI];

                    vertLabels.append(f.size());

                    // If the two normals are colinear, then that means they
                    // point in the same direction
                    if(pos(faceNormal & surfaceNormal))
                    {
                        vtkWriteOps::insert(f + offset, vertLabels);
                    }
                    else
                    {
                        // Need to reverse the face order
                        vtkWriteOps::insert(f.reverseFace() + offset, vertLabels);
                    }
                }
            }
            else
            {
                //Use the original code
                forAll(pp, faceI)
                {
                    const face& f = pp.localFaces()[faceI];

                    vertLabels.append(f.size());
                    vtkWriteOps::insert(f + offset, vertLabels);
                }
            }
            offset += pp.nPoints();
        }
    }

    vtkWriteOps::write(os_, binary_, vertLabels);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Note: The IDs used here are different from the ones written
// by patchWriter::writePatchIDs()
void Foam::patchNFaceSetsWriter::writePatchNFaceSetIDs()
{
    const fvMesh& mesh = vMesh_.mesh();

    DynamicList<floatScalar> fField(nFaces_);

    os_ << "patchID 1 " << nFaces_ << " float" << std::endl;

    forAll(patchIDs_, i)
    {
        label patchI = patchIDs_[i];

        const polyPatch& pp = mesh.boundaryMesh()[patchI];

        if (!isA<emptyPolyPatch>(pp))
        {
            vtkWriteOps::insert(scalarField(pp.size(), i), fField);
        }
    }

    forAll(faceSetIDs_, i)
    {
        label patchI = faceSetIDs_[i];

        if(patchI != -1)
        {
            const primitiveFacePatch& pp = faceSetsAsPrimitives_[i];

            vtkWriteOps::insert(scalarField(pp.size(), patchI), fField);
        }
    }

    vtkWriteOps::write(os_, binary_, fField);
}


// ************************************************************************* //
