/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
 2015-11-10 FSD blueCAPE Lda: Hack for writing faceSets the way we want.
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

#include "manualFaceSetsWriter.H"
#include "vtkWriteFieldOps.H"
#include "faceSet.H"
#include "stringListOps.H"

#if defined( WIN32 ) || defined( WIN64 )
#include "IOobjectList.T.H"
#else
#include "IOobjectList.H"
#endif

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::manualFaceSetsWriter::manualFaceSetsWriter
(
    const vtkMesh& vMesh,
    const bool binary,
    const fileName& fName,
    const PtrList<dictionary>& manuallyTaggedFaceSetsList
)
:
    vMesh_(vMesh),
    binary_(binary),
    fName_(fName),
    faceSetNameIDs_(),
    faceSetCellZoneIDs_(),
    faceSetsAsPrimitives_(),
    os_(fName.c_str(),
        std::ios_base::out|std::ios_base::binary) //a must for Windows!
{
    const fvMesh& mesh = vMesh_.mesh();
    const faceList& meshFaces = mesh.faces();
    const pointField& meshPoints = mesh.points();

    PtrList<const faceSet> faceSets;

    // Read sets
    // Note: Read only the sets we want
    {
        Info<< "Loading faceSets:" << endl;

        faceSetNameIDs_.setSize(manuallyTaggedFaceSetsList.size());
        faceSetCellZoneIDs_.setSize(manuallyTaggedFaceSetsList.size());
        faceSets.setSize(manuallyTaggedFaceSetsList.size());
        faceSetsAsPrimitives_.setSize(manuallyTaggedFaceSetsList.size());

        forAll(manuallyTaggedFaceSetsList, i)
        {
            const dictionary& dict = manuallyTaggedFaceSetsList[i];
            const word setName(dict.lookup("faceSetName"));
            const label nameID(dict.lookupOrDefault("nameID",-1));
            const label cellZoneID(dict.lookupOrDefault("cellZoneID",-1));

            faceSets.set(i, new faceSet(mesh, setName));

            faceSetNameIDs_[i] = nameID;
            faceSetCellZoneIDs_[i] = cellZoneID;

            Info<< "    faceSet " << faceSets[i].name()
                << ", with nameID = " << nameID
                << ", and cellZoneID: " << cellZoneID
                << endl;
        }
    }

    // Write header
    vtkWriteOps::writeHeader(os_, binary_, "patches");
    os_ << "DATASET POLYDATA" << std::endl;

    // Write topology
    nPoints_ = 0;
    nFaces_ = 0;
    label nFaceVerts = 0;

    forAll(faceSets, set)
    {
        const faceSet &fs = faceSets[set];

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

    os_ << "POINTS " << nPoints_ << " float" << std::endl;

    DynamicList<floatScalar> ptField(3*nPoints_);

    forAll(faceSetsAsPrimitives_, i)
    {
        const primitiveFacePatch& pp = faceSetsAsPrimitives_[i];

        vtkWriteOps::insert(pp.localPoints(), ptField);
    }

    vtkWriteOps::write(os_, binary_, ptField);

    os_ << "POLYGONS " << nFaces_ << ' ' << nFaceVerts << std::endl;

    DynamicList<label> vertLabels(nFaceVerts);

    label offset = 0;

    forAll(faceSetsAsPrimitives_, i)
    {
        const primitiveFacePatch& pp = faceSetsAsPrimitives_[i];

        forAll(pp, faceI)
        {
            const face& f = pp.localFaces()[faceI];

            vertLabels.append(f.size());
            vtkWriteOps::insert(f + offset, vertLabels);
        }
        offset += pp.nPoints();
    }

    vtkWriteOps::write(os_, binary_, vertLabels);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Note: The IDs used here are different from the ones written
// by patchWriter::writeNameIDs()
void Foam::manualFaceSetsWriter::writeNameIDs()
{
    DynamicList<label> fField(nFaces_);

    os_ << "nameID 1 " << nFaces_ << " int" << std::endl;

    forAll(faceSetNameIDs_, i)
    {
        label patchI = faceSetNameIDs_[i];

        if(patchI != -1)
        {
            const primitiveFacePatch& pp = faceSetsAsPrimitives_[i];

            vtkWriteOps::insert(labelField(pp.size(), patchI), fField);
        }
    }

    vtkWriteOps::write(os_, binary_, fField);
}


void Foam::manualFaceSetsWriter::writeCellZoneIDs()
{
    DynamicList<label> fField(nFaces_);

    os_ << "cellZoneID 1 " << nFaces_ << " int" << std::endl;

    forAll(faceSetCellZoneIDs_, i)
    {
        label patchI = faceSetCellZoneIDs_[i];

        if(patchI != -1)
        {
            const primitiveFacePatch& pp = faceSetsAsPrimitives_[i];

            vtkWriteOps::insert(labelField(pp.size(), patchI), fField);
        }
    }

    vtkWriteOps::write(os_, binary_, fField);
}


// ************************************************************************* //
