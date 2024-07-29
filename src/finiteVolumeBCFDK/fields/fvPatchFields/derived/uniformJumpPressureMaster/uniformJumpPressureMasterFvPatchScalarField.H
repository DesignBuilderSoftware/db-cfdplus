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
    Foam::uniformJumpPressureMasterFvPatchField

Group
    grpGenericBoundaryConditions

Description
    See the boundary condition uniformJumpPressureFvPatchField for more
    details.

    This class reimplements the jump by adding to the value on the opposite
    patch.

Usage

    \table
        Property        | Description          | Required | Default value
        remotePatchName | name of remote patch | yes      |
        phi             | flux field name      | no       | phi
        specifiedJump   | jump to be applied   | yes      |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            uniformJumpPressureMaster;
        remotePatchName rPatch;
        phi             phi;
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
    uniformJumpPressureMasterFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef uniformJumpPressureMasterFvPatchScalarField_H
#define uniformJumpPressureMasterFvPatchScalarField_H

#include "uniformJumpPressureFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class uniformJumpPressureMasterFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class uniformJumpPressureMasterFvPatchScalarField
:
    public uniformJumpPressureFvPatchField<scalar>
{

protected:

    //- Apply jump
    virtual void applyJump
    (
        const scalar & localAverage,
        const scalar & neighbourAverage,
        const scalar & volumetricFlowRate
    );


public:

    //- Runtime type information
    TypeName("uniformJumpPressureMaster");


    // Constructors

        //- Construct from patch and internal field
        uniformJumpPressureMasterFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        uniformJumpPressureMasterFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given uniformJumpPressureMasterFvPatchScalarField
        //  onto a new patch
        uniformJumpPressureMasterFvPatchScalarField
        (
            const uniformJumpPressureMasterFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        uniformJumpPressureMasterFvPatchScalarField
        (
            const uniformJumpPressureMasterFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new uniformJumpPressureMasterFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        uniformJumpPressureMasterFvPatchScalarField
        (
            const uniformJumpPressureMasterFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar>>
            (
                new uniformJumpPressureMasterFvPatchScalarField(*this, iF)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
