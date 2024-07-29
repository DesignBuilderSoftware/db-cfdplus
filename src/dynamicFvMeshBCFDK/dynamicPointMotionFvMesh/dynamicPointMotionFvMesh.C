/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on a few of OpenFOAM's dynamic meshing classes.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2014-2016 FSD blueCAPE Lda  http://www.bluecape.com.pt
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

#include "dynamicPointMotionFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicPointMotionFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, dynamicPointMotionFvMesh, IOobject);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicPointMotionFvMesh::dynamicPointMotionFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    dynamicMeshCoeffs_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                io.time().constant(),
                *this,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).subDict(typeName + "Coeffs")
    ),
    beforeX_(dynamicMeshCoeffs_.lookup("Xbefore")),
    beforeY_(dynamicMeshCoeffs_.lookup("Ybefore")),
    beforeZ_(dynamicMeshCoeffs_.lookup("Zbefore")),
    afterX_(dynamicMeshCoeffs_.lookup("Xafter")),
    afterY_(dynamicMeshCoeffs_.lookup("Yafter")),
    afterZ_(dynamicMeshCoeffs_.lookup("Zafter")),
    tolerance_(dynamicMeshCoeffs_.lookup<scalar>("tolerance")),
    stationaryPoints_
    (
        IOobject
        (
            "points",
            io.time().constant(),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    )
{
    Info<< "Performing a dynamic mesh displacement from the original axis location: " << endl
        << " Xbefore: " << beforeX_ << endl
        << " Ybefore: " << beforeY_ << endl
        << " Zbefore: " << beforeZ_ << endl
        << "To the new axis locations: " << endl
        << " Xafter: " << afterX_ << endl
        << " Yafter: " << afterY_ << endl
        << " Zafter: " << afterZ_ << endl
        ;

    if(beforeX_.size() != afterX_.size())
    {
        FatalErrorIn("dynamicPointMotionFvMesh(const IOobject& io)")
            << "Xbefore (" << beforeX_.size()
            << ") and Xafter (" << afterX_.size()
            << ") have different lengths."
            << abort(FatalError);
    }
    if(beforeY_.size() != afterY_.size())
    {
        FatalErrorIn("dynamicPointMotionFvMesh(const IOobject& io)")
            << "Ybefore (" << beforeY_.size()
            << ") and Yafter (" << afterY_.size()
            << ") have different lengths."
            << abort(FatalError);
    }
    if(beforeZ_.size() != afterZ_.size())
    {
        FatalErrorIn("dynamicPointMotionFvMesh(const IOobject& io)")
            << "Zbefore (" << beforeZ_.size()
            << ") and Zafter (" << afterZ_.size()
            << ") have different lengths."
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicPointMotionFvMesh::~dynamicPointMotionFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::dynamicPointMotionFvMesh::findSortedIndex
(
    const List<scalar>& l,
    scalar t,
    const label start
)
{
    if (start >= l.size())
    {
        return -1;
    }

    label low = start;
    label high = l.size() - 1;

    while (low <= high)
    {
        label mid = (low + high)/2;

        if (mag(t - l[mid]) < tolerance_)
        {
            return mid;
        }
        else if (t < l[mid])
        {
            high = mid - 1;
        }
        else if (t > l[mid])
        {
            low = mid + 1;
        }
        else
        {
            return mid;
        }
    }

    return -1;
}

Foam::label Foam::dynamicPointMotionFvMesh::findIndex
(
    const List<scalar>& l,
    scalar t,
    const label start
)
{
    label index = -1;

    for (label i = start; i < l.size(); i++)
    {
        if (mag(l[i] - t) < tolerance_)
        {
            index = i;
            break;
        }
    }

    return index;
}

bool Foam::dynamicPointMotionFvMesh::update()
{
    pointField newPoints = stationaryPoints_;
    int countSearchFails(0);

    forAll(newPoints, pointI)
    {
      vector & thePoint = newPoints[pointI];

      if(debug)
      {
        Info << "Before: " << thePoint << endl;
      }

      //lookup the old location
      label xIndex = findSortedIndex(beforeX_, thePoint.x());
      label yIndex = findSortedIndex(beforeY_, thePoint.y());
      label zIndex = findSortedIndex(beforeZ_, thePoint.z());

      //replace with the new location
      if(xIndex >= 0)
      {
          thePoint.x() = afterX_[xIndex];
      }
      else
      {
          countSearchFails++;
      }

      if(yIndex >= 0)
      {
          thePoint.y() = afterY_[yIndex];
      }
      else
      {
          countSearchFails++;
      }

      if(zIndex >= 0)
      {
          thePoint.z() = afterZ_[zIndex];
      }
      else
      {
          countSearchFails++;
      }

      if(debug)
      {
        Info
            << "Indexes found: "
            << xIndex << " "
            << yIndex << " "
            << zIndex << endl;

        Info << "After: " << thePoint << endl << endl;
      }
    }

    fvMesh::movePoints(newPoints);

    Info << "Total number of failed searches: " << countSearchFails << endl;

    //Unused example from dynamicInkJetFvMesh
    /*
    volVectorField& U =
        const_cast<volVectorField&>(lookupObject<volVectorField>("U"));
    U.correctBoundaryConditions();
    */

    return true;
}


// ************************************************************************* //
