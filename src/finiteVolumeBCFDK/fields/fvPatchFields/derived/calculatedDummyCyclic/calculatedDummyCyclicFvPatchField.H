/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2015-2016 FSD blueCAPE Lda  http://www.bluecape.com.pt
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
    Foam::calculatedDummyCyclicFvPatchField

Group
    grpGenericBoundaryConditions

Note
    - Based on OpenFOAM's calculatedFvPatchField. Simply modified to hack
    into the cyclic patch type and it's not meant to be used in fields that
    are part of equations! The objective is to enforce the storage of our
    calculated values.

Description
    This boundary condition is not designed to be evaluated; it is assumed
    that the value is assigned via field assignment, and not via a call to
    e.g. \c updateCoeffs or \c evaluate.

    \heading Patch usage

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            calculatedDummyCyclic;
        patchType       cyclic;
        value           uniform (0 0 0);    // optional value entry
    }
    \endverbatim

SourceFiles
    calculatedDummyCyclicFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef calculatedDummyCyclicFvPatchField_H
#define calculatedDummyCyclicFvPatchField_H

#include "calculatedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class calculatedDummyCyclicFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class calculatedDummyCyclicFvPatchField
:
    public calculatedFvPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("calculatedDummyCyclic");


    // Constructors

        //- Construct from patch and internal field
        calculatedDummyCyclicFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        calculatedDummyCyclicFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        calculatedDummyCyclicFvPatchField
        (
            const calculatedDummyCyclicFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        calculatedDummyCyclicFvPatchField
        (
            const calculatedDummyCyclicFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new calculatedDummyCyclicFvPatchField<Type>(*this)
            );
        }

        //- Copy constructor setting internal field reference
        calculatedDummyCyclicFvPatchField
        (
            const calculatedDummyCyclicFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new calculatedDummyCyclicFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "calculatedDummyCyclicFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
