/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::uniformJumpPressureFvPatchField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a uniform fixed jump value condition,
    which works in function of the average going flow rate accross the local
    patch. The jump is subtracted from the average value at the remote patch.

    This boundary condition was designed for the pressure field, to apply a
    pressure jump on the upstream side, however this can be applied to the
    majority of the fields that may need a jump in function of the flow rate.

    The result of the evaluated jump is expected to be in the same units as the
    field this boundary is defined, it does not convert based on density.

Usage

    \table
        Property        | Description          | Required | Default value
        remotePatchName | name of remote patch | yes      |
        phi             | flux field name      | no       | phi
        rho             | density field name   | no       | rho
        specifiedJump   | jump to be applied   | yes      |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            uniformJumpPressure;
        remotePatchName rPatch;
        phi             phi;
        rho             none; //for incompressible flow
        specifiedJump   constant 1;
        relaxation      0.5;
        value           uniform 0; // place holder
    }
    \endverbatim

Note
    The specifiedJump entry is a Function1 type, able to describe time
    varying functions. The example above gives the usage for supplying a
    constant value.

    This class was based on OpenFOAM's uniformFixedValueFvPatchField and
    outletMappedUniformInletFvPatchField.

See also
    Foam::Function1s
    Foam::fixedValueFvPatchField

SourceFiles
    uniformJumpPressureFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef uniformJumpPressureFvPatchField_H
#define uniformJumpPressureFvPatchField_H

#include "fixedValueFvPatchFields.H"

#if defined( WIN32 ) || defined( WIN64 )
#include "Function1.T.H"
#else
#include "Function1.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class uniformJumpPressureFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class uniformJumpPressureFvPatchField
:
    public fixedValueFvPatchField<Type>
{

protected:

    // Protected data

        //- Name of the remote patch to be mapped
        word remotePatchName_;

        //- Name of the flux transporting the field
        word phiName_;

        //- Name of the density field used to normalize the mass flux
        word rhoName_;

        //- Variable jump specified by the user
        autoPtr<Function1<Type>> userSpecifiedJump_;

        //- Relaxation factor
        scalar relaxation_;


    //- Apply jump
    virtual void applyJump
    (
        const Type & localAverage,
        const Type & neighbourAverage,
        const scalar & volumetricFlowRate
    );


public:

    //- Runtime type information
    TypeName("uniformJumpPressure");


    // Constructors

        //- Construct from patch and internal field
        uniformJumpPressureFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch and internal field and patch field
        uniformJumpPressureFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const Field<Type>& fld
        );

        //- Construct from patch, internal field and dictionary
        uniformJumpPressureFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given uniformJumpPressureFvPatchField
        //  onto a new patch
        uniformJumpPressureFvPatchField
        (
            const uniformJumpPressureFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        uniformJumpPressureFvPatchField
        (
            const uniformJumpPressureFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new uniformJumpPressureFvPatchField<Type>(*this)
            );
        }

        //- Copy constructor setting internal field reference
        uniformJumpPressureFvPatchField
        (
            const uniformJumpPressureFvPatchField<Type>&,
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
                new uniformJumpPressureFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Name of the remote patch to be mapped
            const word& remotePatchName() const
            {
                return remotePatchName_;
            }

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "uniformJumpPressureFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
