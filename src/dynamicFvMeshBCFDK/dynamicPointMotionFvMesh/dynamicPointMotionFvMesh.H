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

Class
    Foam::dynamicPointMotionFvMesh

Description
    Mesh motion specifically for the change of all point locations from one
    general position to another, using uniform before and after locations
    along the major axis.

    Warning: This class is only meant to perform mesh motion and will not
    do any specific field corrections!

    Example:

    \verbatim
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        object      dynamicMeshDict;
    }

    dynamicFvMeshLibs  ("libdynamicFvMeshBCFDK.so");

    dynamicFvMesh     dynamicPointMotionFvMesh;

    motionSolverLibs   ("libfvMotionSolvers.so");

    dynamicPointMotionFvMeshCoeffs
    {
      Xbefore (0.0 0.5 1.0);
      Ybefore (0.0 0.2 0.4);
      Zbefore (0.0 1.0 2.0);

      Xafter (0.0 0.35 1.0);
      Yafter (0.0 0.3 0.4);
      Zafter (0.0 1.5 2.0);
    }
    \endverbatim


SourceFiles
    dynamicPointMotionFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicPointMotionFvMesh_H
#define dynamicPointMotionFvMesh_H

#include "dynamicFvMesh.H"
#include "dictionary.H"
#include "pointIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dynamicPointMotionFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicPointMotionFvMesh
:
    public dynamicFvMesh
{
    // Private Data

        dictionary dynamicMeshCoeffs_;

        List<scalar> beforeX_;
        List<scalar> beforeY_;
        List<scalar> beforeZ_;

        List<scalar> afterX_;
        List<scalar> afterY_;
        List<scalar> afterZ_;

        scalar tolerance_;

        pointIOField stationaryPoints_;


    // Private Member Functions

        label findSortedIndex
        (
            const List<scalar>& l,
            scalar t,
            const label start=0
        );

        label findIndex
        (
            const List<scalar>& l,
            scalar t,
            const label start=0
        );


public:

    //- Runtime type information
    TypeName("dynamicPointMotionFvMesh");


    // Constructors

        //- Construct from IOobject
        dynamicPointMotionFvMesh(const IOobject& io);

        //- Disallow default bitwise copy construction
        dynamicPointMotionFvMesh(const dynamicPointMotionFvMesh&) = delete;

    //- Destructor
    ~dynamicPointMotionFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const dynamicPointMotionFvMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
