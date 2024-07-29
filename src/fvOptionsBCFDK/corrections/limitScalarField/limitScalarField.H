/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's limitVelocity.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2021 FSD blueCAPE Lda  http://www.bluecape.com.pt
-------------------------------------------------------------------------------
License
    This file is a derivative of OpenFOAM.

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
    Foam::fv::limitScalarField

Description
    Limits the maximum velocity magnitude to the specified \c max value.

Usage
    Example usage:
    \verbatim
    limitU
    {
        type            limitScalarField;
        active          yes;

        selectionMode   all;
        max             100;
    }
    \endverbatim

SourceFiles
    limitScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef limitScalarField_H
#define limitScalarField_H

#include "cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class limitScalarField Declaration
\*---------------------------------------------------------------------------*/

class limitScalarField
:
    public cellSetOption
{
protected:

    // Protected data

        //- Limited field name
        word targetFieldName_;

        //- Reference field name to get maximum from
        word refFieldName_;

        //- Maximum scale up value
        scalar scaleMax_;


public:

    //- Runtime type information
    TypeName("limitScalarField");


    // Constructors

        //- Construct from components
        limitScalarField
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        limitScalarField(const limitScalarField&) = delete;


    //- Destructor
    virtual ~limitScalarField()
    {}


    // Member Functions

        //- Read dictionary
        virtual bool read(const dictionary& dict);

        //- Correct the target field
        virtual void correct(volScalarField& field);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const limitScalarField&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
