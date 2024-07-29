/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2014-2018 FSD blueCAPE Lda: Added the method writeCellZoneIDs()
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

#include "internalWriterBCFDK.H"
#include "vtkWriteFieldOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::internalWriterBCFDK::internalWriterBCFDK
(
    const vtkMesh& vMesh,
    const bool binary,
    const fileName& fName
)
:
    internalWriter(vMesh, binary, fName),
    vMesh_(vMesh),
    binary_(binary)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::internalWriterBCFDK::writeCellZoneIDs()
{
    const fvMesh& mesh = vMesh_.mesh();
    const vtkTopo& topo = vMesh_.topo();
    const labelList& vtkCellTypes = topo.cellTypes();
    const labelList& superCells = topo.superCells();
    const cellZoneMesh &czMesh = mesh.cellZones();
    labelList czIndexToID = labelList(czMesh.size(), -1);

    //give a list of the existing cell zones
    {
        Info<< "Exporting the following IDs for the cellZoneID field:" << endl;
        Info<< "    internalMesh:   -1" << endl;

        label counter(-1);
        forAll(czMesh, i)
        {
            if(!czMesh[i].empty())
            {
              counter++;
              czIndexToID[i] = counter;

              Info<< "    " << czMesh[i].name() << ":   " << counter << endl;
            }
        }
    }

    // Cell ids first
    os() << "cellZoneID 1 " << vtkCellTypes.size() << " int" << std::endl;

    labelList cellId(vtkCellTypes.size());
    label labelI = 0;

    if (vMesh_.useSubMesh())
    {
        //WARNING: This part of the code hasn't been tested yet!!
        const labelList& cMap = vMesh_.subsetter().cellMap();

        forAll(mesh.cells(), cellI)
        {
            label zoneID = czMesh.whichZone(cMap[cellI]);
            if(zoneID>=0)
            {
              zoneID = czIndexToID[zoneID];
            }
            cellId[labelI++] = zoneID;
        }
        forAll(superCells, superCellI)
        {
            label origCellI = cMap[superCells[superCellI]];

            label zoneID = czMesh.whichZone(origCellI);
            if(zoneID>=0)
            {
              zoneID = czIndexToID[zoneID];
            }
            cellId[labelI++] = zoneID;
        }
    }
    else
    {
        forAll(mesh.cells(), cellI)
        {
            label zoneID = czMesh.whichZone(cellI);
            if(zoneID>=0)
            {
              zoneID = czIndexToID[zoneID];
            }
            cellId[labelI++] = zoneID;
        }
        forAll(superCells, superCellI)
        {
            label origCellI = superCells[superCellI];

            label zoneID = czMesh.whichZone(origCellI);
            if(zoneID>=0)
            {
              zoneID = czIndexToID[zoneID];
            }
            cellId[labelI++] = zoneID;
        }
    }

    vtkWriteOps::write(os(), binary_, cellId);
}



// ************************************************************************* //
