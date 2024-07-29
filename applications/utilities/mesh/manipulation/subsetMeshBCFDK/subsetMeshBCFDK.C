/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's subsetMesh.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2015-2018 FSD blueCAPE Lda  http://www.bluecape.com.pt
-------------------------------------------------------------------------------
License
    This file is a derivative work from OpenFOAM.

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
    subsetMeshBCFDK

Description
    Selects a section of mesh based on a cellSet.

    The utility sub-sets the mesh to choose only a part of interest. Check
    the setSet/cellSet/topoSet utilities to see how to select cells based on
    various shapes.

    The mesh will subset all points, faces and cells needed to make a sub-mesh
    but will not preserve attached boundary types.

    Modified subsetMesh from OpenFOAM 5.x, to also support remapping of
    faceSets. cellSets and pointSets are not supported yet, since we don't
    need it so far.

\*---------------------------------------------------------------------------*/

#include "fvMeshSubset.H"
#include "argList.H"
#include "cellSet.H"
#include "faceSet.H"
#include "labelToFace.H"
#include "volFields.H"

#if defined(WIN32) || defined(WIN64)
#include "IOobjectList.T.H"
#else
#include "IOobjectList.H"
#endif

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


template<class Type>
void subsetVolFields
(
    const fvMeshSubset& subsetter,
    const wordList& fieldNames,
    PtrList<GeometricField<Type, fvPatchField, volMesh>>& subFields
)
{
    const fvMesh& baseMesh = subsetter.baseMesh();

    forAll(fieldNames, i)
    {
        const word& fieldName = fieldNames[i];

        Info<< "Subsetting field " << fieldName << endl;

        GeometricField<Type, fvPatchField, volMesh> fld
        (
            IOobject
            (
                fieldName,
                baseMesh.time().timeName(),
                baseMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            baseMesh
        );

        subFields.set(i, subsetter.interpolate(fld));
    }
}


template<class Type>
void subsetSurfaceFields
(
    const fvMeshSubset& subsetter,
    const wordList& fieldNames,
    PtrList<GeometricField<Type, fvsPatchField, surfaceMesh>>& subFields
)
{
    const fvMesh& baseMesh = subsetter.baseMesh();

    forAll(fieldNames, i)
    {
        const word& fieldName = fieldNames[i];

        Info<< "Subsetting field " << fieldName << endl;

        GeometricField<Type, fvsPatchField, surfaceMesh> fld
        (
            IOobject
            (
                fieldName,
                baseMesh.time().timeName(),
                baseMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            baseMesh
        );

        subFields.set(i, subsetter.interpolate(fld));
    }
}


template<class Type>
void subsetPointFields
(
    const fvMeshSubset& subsetter,
    const pointMesh& pMesh,
    const wordList& fieldNames,
    PtrList<GeometricField<Type, pointPatchField, pointMesh>>& subFields
)
{
    const fvMesh& baseMesh = subsetter.baseMesh();

    forAll(fieldNames, i)
    {
        const word& fieldName = fieldNames[i];

        Info<< "Subsetting field " << fieldName << endl;

        GeometricField<Type, pointPatchField, pointMesh> fld
        (
            IOobject
            (
                fieldName,
                baseMesh.time().timeName(),
                baseMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            pMesh
        );

        subFields.set(i, subsetter.interpolate(fld));
    }
}


template<class Type>
void subsetDimensionedFields
(
    const fvMeshSubset& subsetter,
    const wordList& fieldNames,
    PtrList<DimensionedField<Type, volMesh>>& subFields
)
{
    const fvMesh& baseMesh = subsetter.baseMesh();

    forAll(fieldNames, i)
    {
        const word& fieldName = fieldNames[i];

        Info<< "Subsetting field " << fieldName << endl;

        DimensionedField<Type, volMesh> fld
        (
            IOobject
            (
                fieldName,
                baseMesh.time().timeName(),
                baseMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            baseMesh
        );

        subFields.set(i, subsetter.interpolate(fld));
    }
}



void populateReverseFaceMap
(
    const fvMesh& mesh,
    const labelList& faceMap,
    labelList &reverseFaceMap
)
{
    reverseFaceMap.setSize(mesh.nFaces(), -1);

    forAll(faceMap, i)
    {
        if(faceMap[i]>=0)
        {
            reverseFaceMap[faceMap[i]] = i;
        }
    }
}


void remapFaceSet
(
    const fvMesh& subMesh,
    faceSet& faceSetToChange,
    const labelList& reverseFaceMap
)
{
    labelList newFaceIdList;
    label newIDCounter = 0;

    newFaceIdList.setSize(faceSetToChange.size());

    forAllConstIter(labelHashSet, faceSetToChange, iter)
    {
        const label newFaceI = reverseFaceMap[iter.key()];

        if (newFaceI >= 0)
        {
            newFaceIdList[newIDCounter]=newFaceI;
            newIDCounter++;
        }
    }

    //reset size
    newFaceIdList.setSize(newIDCounter);

    labelToFace source
    (
        subMesh,
        newFaceIdList
    );

    faceSetToChange.clearStorage();
    source.applyToSet(topoSetSource::NEW, faceSetToChange);
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "select a mesh subset based on a cellSet"
    );

    #include "addOverwriteOption.H"
    #include "addRegionOption.H"
    argList::validArgs.append("cellSet");
    argList::addOption
    (
        "patch",
        "name",
        "add exposed internal faces to specified patch instead of to "
        "'oldInternalFaces'"
    );
    argList::addOption
    (
        "resultTime",
        "time",
        "specify a time for the resulting mesh"
    );
    argList::addBoolOption
    (
        "remapFaceSets",
        "also remap all faceSets from the old mesh to the new mesh"
    );
    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();

    Foam::word meshRegionName = polyMesh::defaultRegion;
    args.optionReadIfPresent("region", meshRegionName);

    #include "createNamedMesh.H"


    const word setName = args[1];

    word meshInstance = mesh.pointsInstance();
    word fieldsInstance = runTime.timeName();

    const bool overwrite = args.optionFound("overwrite");
    const bool specifiedInstance = args.optionReadIfPresent
    (
        "resultTime",
        fieldsInstance
    );
    if (specifiedInstance)
    {
        // Set both mesh and field to this time
        meshInstance = fieldsInstance;
    }

    const bool remapFaceSets = args.optionFound("remapFaceSets");


    Info<< "Reading cell set from " << setName << endl << endl;

    // Create mesh subsetting engine
    fvMeshSubset subsetter(mesh);

    label patchi = -1;

    if (args.optionFound("patch"))
    {
        const word patchName = args["patch"];

        patchi = mesh.boundaryMesh().findPatchID(patchName);

        if (patchi == -1)
        {
            FatalErrorInFunction
                << nl << "Valid patches are " << mesh.boundaryMesh().names()
                << exit(FatalError);
        }

        Info<< "Adding exposed internal faces to patch " << patchName << endl
            << endl;
    }
    else
    {
        Info<< "Adding exposed internal faces to a patch called"
            << " \"oldInternalFaces\" (created if necessary)" << endl
            << endl;
    }


    cellSet currentSet(mesh, setName);

    subsetter.setLargeCellSubset(currentSet, patchi, true);

    IOobjectList objects(mesh, runTime.timeName());


    // Read vol fields and subset
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    wordList scalarNames(objects.names(volScalarField::typeName));
    PtrList<volScalarField> scalarFlds(scalarNames.size());
    subsetVolFields(subsetter, scalarNames, scalarFlds);

    wordList vectorNames(objects.names(volVectorField::typeName));
    PtrList<volVectorField> vectorFlds(vectorNames.size());
    subsetVolFields(subsetter, vectorNames, vectorFlds);

    wordList sphericalTensorNames
    (
        objects.names(volSphericalTensorField::typeName)
    );
    PtrList<volSphericalTensorField> sphericalTensorFlds
    (
        sphericalTensorNames.size()
    );
    subsetVolFields(subsetter, sphericalTensorNames, sphericalTensorFlds);

    wordList symmTensorNames(objects.names(volSymmTensorField::typeName));
    PtrList<volSymmTensorField> symmTensorFlds(symmTensorNames.size());
    subsetVolFields(subsetter, symmTensorNames, symmTensorFlds);

    wordList tensorNames(objects.names(volTensorField::typeName));
    PtrList<volTensorField> tensorFlds(tensorNames.size());
    subsetVolFields(subsetter, tensorNames, tensorFlds);


    // Read surface fields and subset
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    wordList surfScalarNames(objects.names(surfaceScalarField::typeName));
    PtrList<surfaceScalarField> surfScalarFlds(surfScalarNames.size());
    subsetSurfaceFields(subsetter, surfScalarNames, surfScalarFlds);

    wordList surfVectorNames(objects.names(surfaceVectorField::typeName));
    PtrList<surfaceVectorField> surfVectorFlds(surfVectorNames.size());
    subsetSurfaceFields(subsetter, surfVectorNames, surfVectorFlds);

    wordList surfSphericalTensorNames
    (
        objects.names(surfaceSphericalTensorField::typeName)
    );
    PtrList<surfaceSphericalTensorField> surfSphericalTensorFlds
    (
        surfSphericalTensorNames.size()
    );
    subsetSurfaceFields
    (
        subsetter,
        surfSphericalTensorNames,
        surfSphericalTensorFlds
    );

    wordList surfSymmTensorNames
    (
        objects.names(surfaceSymmTensorField::typeName)
    );
    PtrList<surfaceSymmTensorField> surfSymmTensorFlds
    (
        surfSymmTensorNames.size()
    );
    subsetSurfaceFields(subsetter, surfSymmTensorNames, surfSymmTensorFlds);

    wordList surfTensorNames(objects.names(surfaceTensorField::typeName));
    PtrList<surfaceTensorField> surfTensorFlds(surfTensorNames.size());
    subsetSurfaceFields(subsetter, surfTensorNames, surfTensorFlds);


    // Read point fields and subset
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const pointMesh& pMesh = pointMesh::New(mesh);

    wordList pointScalarNames(objects.names(pointScalarField::typeName));
    PtrList<pointScalarField> pointScalarFlds(pointScalarNames.size());
    subsetPointFields(subsetter, pMesh, pointScalarNames, pointScalarFlds);

    wordList pointVectorNames(objects.names(pointVectorField::typeName));
    PtrList<pointVectorField> pointVectorFlds(pointVectorNames.size());
    subsetPointFields(subsetter, pMesh, pointVectorNames, pointVectorFlds);

    wordList pointSphericalTensorNames
    (
        objects.names(pointSphericalTensorField::typeName)
    );
    PtrList<pointSphericalTensorField> pointSphericalTensorFlds
    (
        pointSphericalTensorNames.size()
    );
    subsetPointFields
    (
        subsetter,
        pMesh,
        pointSphericalTensorNames,
        pointSphericalTensorFlds
    );

    wordList pointSymmTensorNames
    (
        objects.names(pointSymmTensorField::typeName)
    );
    PtrList<pointSymmTensorField> pointSymmTensorFlds
    (
        pointSymmTensorNames.size()
    );
    subsetPointFields
    (
        subsetter,
        pMesh,
        pointSymmTensorNames,
        pointSymmTensorFlds
    );

    wordList pointTensorNames(objects.names(pointTensorField::typeName));
    PtrList<pointTensorField> pointTensorFlds(pointTensorNames.size());
    subsetPointFields(subsetter, pMesh, pointTensorNames, pointTensorFlds);


    // Read dimensioned fields and subset
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    typedef volScalarField::Internal dimScalType;
    wordList scalarDimNames(objects.names(dimScalType::typeName));
    PtrList<dimScalType> scalarDimFlds(scalarDimNames.size());
    subsetDimensionedFields(subsetter, scalarDimNames, scalarDimFlds);

    typedef volVectorField::Internal dimVecType;
    wordList vectorDimNames(objects.names(dimVecType::typeName));
    PtrList<dimVecType> vectorDimFlds(vectorDimNames.size());
    subsetDimensionedFields(subsetter, vectorDimNames, vectorDimFlds);

    typedef volSphericalTensorField::Internal dimSphereType;
    wordList sphericalTensorDimNames(objects.names(dimSphereType::typeName));
    PtrList<dimSphereType> sphericalTensorDimFlds
    (
        sphericalTensorDimNames.size()
    );
    subsetDimensionedFields
    (
        subsetter,
        sphericalTensorDimNames,
        sphericalTensorDimFlds
    );

    typedef volSymmTensorField::Internal dimSymmTensorType;
    wordList symmTensorDimNames(objects.names(dimSymmTensorType::typeName));
    PtrList<dimSymmTensorType> symmTensorDimFlds(symmTensorDimNames.size());
    subsetDimensionedFields(subsetter, symmTensorDimNames, symmTensorDimFlds);

    typedef volTensorField::Internal dimTensorType;
    wordList tensorDimNames(objects.names(dimTensorType::typeName));
    PtrList<dimTensorType> tensorDimFlds(tensorDimNames.size());
    subsetDimensionedFields(subsetter, tensorDimNames, tensorDimFlds);


    // Write mesh and fields to new time
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (overwrite || specifiedInstance)
    {
        runTime.setTime(instant(fieldsInstance), 0);
        subsetter.subMesh().setInstance(meshInstance);
    }
    else
    {
        runTime++;
    }

    Info<< "Writing subsetted mesh and fields to time " << runTime.timeName()
        << endl;
    subsetter.subMesh().write();


    // Subsetting adds 'subset' prefix. Rename field to be like original.
    forAll(scalarFlds, i)
    {
        scalarFlds[i].rename(scalarNames[i]);
        scalarFlds[i].write();
    }
    forAll(vectorFlds, i)
    {
        vectorFlds[i].rename(vectorNames[i]);
        vectorFlds[i].write();
    }
    forAll(sphericalTensorFlds, i)
    {
        sphericalTensorFlds[i].rename(sphericalTensorNames[i]);
        sphericalTensorFlds[i].write();
    }
    forAll(symmTensorFlds, i)
    {
        symmTensorFlds[i].rename(symmTensorNames[i]);
        symmTensorFlds[i].write();
    }
    forAll(tensorFlds, i)
    {
        tensorFlds[i].rename(tensorNames[i]);
        tensorFlds[i].write();
    }

    // Surface ones.
    forAll(surfScalarFlds, i)
    {
        surfScalarFlds[i].rename(surfScalarNames[i]);
        surfScalarFlds[i].write();
    }
    forAll(surfVectorFlds, i)
    {
        surfVectorFlds[i].rename(surfVectorNames[i]);
        surfVectorFlds[i].write();
    }
    forAll(surfSphericalTensorFlds, i)
    {
        surfSphericalTensorFlds[i].rename(surfSphericalTensorNames[i]);
        surfSphericalTensorFlds[i].write();
    }
    forAll(surfSymmTensorFlds, i)
    {
        surfSymmTensorFlds[i].rename(surfSymmTensorNames[i]);
        surfSymmTensorFlds[i].write();
    }
    forAll(surfTensorNames, i)
    {
        surfTensorFlds[i].rename(surfTensorNames[i]);
        surfTensorFlds[i].write();
    }

    // Point ones
    forAll(pointScalarFlds, i)
    {
        pointScalarFlds[i].rename(pointScalarNames[i]);
        pointScalarFlds[i].write();
    }
    forAll(pointVectorFlds, i)
    {
        pointVectorFlds[i].rename(pointVectorNames[i]);
        pointVectorFlds[i].write();
    }
    forAll(pointSphericalTensorFlds, i)
    {
        pointSphericalTensorFlds[i].rename(pointSphericalTensorNames[i]);
        pointSphericalTensorFlds[i].write();
    }
    forAll(pointSymmTensorFlds, i)
    {
        pointSymmTensorFlds[i].rename(pointSymmTensorNames[i]);
        pointSymmTensorFlds[i].write();
    }
    forAll(pointTensorNames, i)
    {
        pointTensorFlds[i].rename(pointTensorNames[i]);
        pointTensorFlds[i].write();
    }

    // DimensionedFields
    forAll(scalarDimFlds, i)
    {
        scalarDimFlds[i].rename(scalarDimNames[i]);
        scalarDimFlds[i].write();
    }
    forAll(vectorDimFlds, i)
    {
        vectorDimFlds[i].rename(vectorDimNames[i]);
        vectorDimFlds[i].write();
    }
    forAll(sphericalTensorDimFlds, i)
    {
        sphericalTensorDimFlds[i].rename(sphericalTensorDimNames[i]);
        sphericalTensorDimFlds[i].write();
    }
    forAll(symmTensorDimFlds, i)
    {
        symmTensorDimFlds[i].rename(symmTensorDimNames[i]);
        symmTensorDimFlds[i].write();
    }
    forAll(tensorDimFlds, i)
    {
        tensorDimFlds[i].rename(tensorDimNames[i]);
        tensorDimFlds[i].write();
    }

    if (remapFaceSets)
    {
        // Read faceSets
        PtrList<faceSet> faceSets;
        labelList reverseFaceMap;

        Info<< "Generating reverse face map..." << endl;

        populateReverseFaceMap(mesh, subsetter.faceMap(), reverseFaceMap);

        Info<< "Loading faceSets for remapping..." << endl;

        // Read sets
        IOobjectList objects(mesh, mesh.facesInstance(), "polyMesh/sets");
        {
            IOobjectList fSets(objects.lookupClass(faceSet::typeName));
            if (fSets.size())
            {
                label setCounter = 0;
                faceSets.setSize(fSets.size());

                Info<< "Reading faceSets:" << endl;
                forAllConstIter(IOobjectList, fSets, iter)
                {
                    faceSets.set(setCounter, new faceSet(*iter()));
                    Info<< "    " << faceSets[setCounter].name() << endl;

                    setCounter++;
                }
            }
        }

        // Remap faceSets to the new mesh structure
        forAll(faceSets, i)
        {
            remapFaceSet(subsetter.subMesh(), faceSets[i], reverseFaceMap);
            faceSets[i].instance() = subsetter.subMesh().facesInstance();

            //write with stricter options, because the default messes things a
            //bit, and mistags the header with binary when it's ASCII.
            faceSets[i].writeObject(
                IOstream::ASCII,
                mesh.time().writeVersion(),
                IOstream::UNCOMPRESSED,
                true
                );
            //faceSets[i].write();
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
