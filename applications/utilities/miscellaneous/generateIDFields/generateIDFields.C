/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on various OpenFOAM utilities.
    Developed for blueCFD(R)-Kernel
    Copyright (C) 2015-2018 FSD blueCAPE Lda  http://www.bluecape.com.pt
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
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Application
    generateIDFields

Description
    Generates the following identification fields:
      - cellID for cells
      - sideID for patches
      - Normal for faces on the surface mesh.

    Depends on the generateIDFieldsDict dictionary file.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"

#if defined( WIN32 ) || defined( WIN64 )
#include "Time.T.H"
#else
#include "Time.H"
#endif

#include "fvMesh.H"
#include "vectorIOField.H"
#include "volFields.H"
#include "cyclicFvPatch.H"

//Hack so that we don't need to link directly to our finiteVolumeBCFDK
#define CALCULATED_DUMMY_CYCLIC_TYPENAME "calculatedDummyCyclic"

int main(int argc, char *argv[])
{

    timeSelector::addOptions();
#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createNamedMesh.H"

   /**
     * IODictionary: generateIDFieldsDict
     * Purpose: For now, provides the list of patches that are tagged
     * as sideID = 1.
     */
    IOdictionary generateIDFieldsDict
    (
        IOobject
        (
            "generateIDFieldsDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    /**
     * Word list for the names of the patches to be tagged as sideID = 1
     */
    const wordList patchSideBNames(
        generateIDFieldsDict.lookup("sideBPatchNames")
        );

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        // Check for new mesh
        mesh.readUpdate();

        //------------------------------

        Info<< "Writing cellID" << endl;

        {
            volScalarField cellID
            (
                IOobject
                (
                    "cellID",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("cellID", dimless, scalar(0))
            );

            {
                forAll(cellID, itemI)
                {
                  cellID[itemI] = scalar(itemI);
                }

                cellID.write();
            }
        }

        //------------------------------

        const fvBoundaryMesh& bm = mesh.boundary();

        wordList allBoundaryTypesAsCalculated(bm.size());
        wordList actualBoundaryTypes(bm.size());

        forAll(bm, patchI)
        {
            if (isA<cyclicFvPatch>(bm[patchI]))
            {
                actualBoundaryTypes[patchI] =
                    cyclicFvPatch::typeName;

                allBoundaryTypesAsCalculated[patchI] =
                    CALCULATED_DUMMY_CYCLIC_TYPENAME;
            }
            //Note: There is no need to specifically check all other types of
            //patches, because OpenFOAM will only apply the patch types that
            //are compatible. At least in OpenFOAM 2.2.x
            else
            {
                actualBoundaryTypes[patchI] =
                    calculatedFvPatchField<scalar>::typeName;

                allBoundaryTypesAsCalculated[patchI] =
                    calculatedFvPatchField<scalar>::typeName;
            }

        }

        //------------------------------

        Info<< "Writing SurfaceNormal" << endl;

        {
            volVectorField SurfaceNormal
            (
                IOobject
                (
                    "SurfaceNormal",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedVector("SurfaceNormal", dimless, vector(0,0,0)),
                allBoundaryTypesAsCalculated,
                actualBoundaryTypes
            );

            {
                const fvPatchList& patches = mesh.boundary();

                forAll(patches, patchI)
                {
                    // We want the normal to point into the domain,
                    // hence the minus
                    SurfaceNormal.boundaryFieldRef()[patchI] =
                        -patches[patchI].nf();
                }

                SurfaceNormal.write();
            }
        }

        //------------------------------

        Info<< "Writing sideID" << endl;

        {
            volScalarField sideID
            (
                IOobject
                (
                    "sideID",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("sideID", dimless, scalar(0)),
                allBoundaryTypesAsCalculated,
                actualBoundaryTypes
            );

            {
                forAll(patchSideBNames, bnameI)
                {
                    label patchI = mesh.boundaryMesh().
                        findPatchID(patchSideBNames[bnameI]);

                    if(patchI>=0)
                    {
                        fvPatchField<scalar>& pfld =
                            sideID.boundaryFieldRef()[patchI];

                        forAll(pfld, faceI)
                        {
                            pfld[faceI] = scalar(1.0); //side B = 1
                        }
                    }
                }

                sideID.write();
            }
        }
    }

    Info << "End\n" << endl;

    return 0;
}

