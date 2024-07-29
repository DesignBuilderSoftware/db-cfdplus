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

Class
    Foam::internalWriterBCFDK

Description
    Inherits from internalWriter to add the new method.

SourceFiles
    internalWriterBCFDK.C
    internalWriterBCFDKTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef internalWriterBCFDK_H
#define internalWriterBCFDK_H

#include "OFstream.H"
#include "volFields.H"
#include "pointFields.H"
#include "vtkMesh.H"
#include "internalWriter.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class internalWriterBCFDK Declaration
\*---------------------------------------------------------------------------*/

class internalWriterBCFDK : public internalWriter
{
    //We need to keep our own copy
    const vtkMesh& vMesh_;

    const bool binary_;

public:

    // Constructors

        //- Construct from components
        internalWriterBCFDK
        (
            const vtkMesh&,
            const bool binary,
            const fileName&
        );


    // Member Functions

        //- Write cellZoneIDs
        void writeCellZoneIDs();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
