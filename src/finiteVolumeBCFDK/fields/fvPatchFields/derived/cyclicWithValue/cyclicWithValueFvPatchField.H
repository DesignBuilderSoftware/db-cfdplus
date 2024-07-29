/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Modified and developed for blueCFD(R)-Kernel
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
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::cyclicWithValueFvPatchField

Group
    grpCoupledBoundaryConditions

Description
    This boundary condition is designed to work in everything as a cyclic
    boundary condition, with the exception it also writes the value.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            cyclicWithValue;
        patchType       cyclic;
        value           uniform (0 0 0);    // optional value entry
    }
    \endverbatim

SourceFiles
    cyclicWithValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicWithValueFvPatchField_H
#define cyclicWithValueFvPatchField_H

#include "cyclicFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class cyclicWithValueFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class cyclicWithValueFvPatchField
:
    public cyclicFvPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("cyclicWithValue");


    // Constructors

        //- Construct from patch and internal field
        cyclicWithValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        cyclicWithValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        cyclicWithValueFvPatchField
        (
            const cyclicWithValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        cyclicWithValueFvPatchField
        (
            const cyclicWithValueFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new cyclicWithValueFvPatchField<Type>(*this)
            );
        }

        //- Copy constructor setting internal field reference
        cyclicWithValueFvPatchField
        (
            const cyclicWithValueFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new cyclicWithValueFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Return the interface type
            virtual const word& interfaceFieldType() const
            {
                return cyclicFvPatchField<Type>::type();
            }

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "cyclicWithValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
