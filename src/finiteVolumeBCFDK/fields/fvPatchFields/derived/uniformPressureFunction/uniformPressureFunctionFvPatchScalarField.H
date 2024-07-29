/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2015-2021 FSD blueCAPE Lda  http://www.bluecape.com.pt
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
    Foam::uniformPressureFunctionFvPatchScalarField

Group
    grpInletBoundaryConditions grpOutletBoundaryConditions

Description
    This boundary condition can be applied to assign either a pressure inlet
    or outlet total pressure condition for a ventilator.

Usage

    \table
        Property     | Description             | Required    | Default value
        function     | flow rate to pressure function to be applied   | yes  |
        direction    | direction of flow through ventilator [in/out] | yes |
        p0           | environmental total pressure | yes    |
        minFlowRate  | minimum flow-rate allowed (>=0) | no | not limited
        maxFlowRate  | maximum flow-rate allowed (>0) | no | not limited
        relaxation   | factor applied to the result ]0,1] | no | 1.0
    \endtable

    Other properties are defined in dynamicPressureFvPatchScalarField.

    Example of the boundary condition specification:
    \verbatim
    inlet
    {
        type            uniformPressureFunction;
        function        constant 1;
        phi             phi;
        direction       in;
        p0              uniform 0;
        value           uniform 0;
    }

    outlet
    {
        type            uniformPressureFunction;
        function        constant 1;
        direction       out;
        p0              uniform 0;
        value           uniform 0;
    }
    \endverbatim

Note
    The function entry is a Function1 type, able to describe time
    varying functions. The example above gives the usage for supplying a
    constant value.

    This class was based on OpenFOAM's fanPressureFvPatchField and based on
    blueCFD-Kernel's class uniformJumpPressureFvPatchField.

See also
    Foam::dynamicPressureFvPatchScalarField
    Foam::uniformJumpPressureFvPatchField
    Foam::Function1

SourceFiles
   uniformPressureFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef uniformPressureFunctionFvPatchScalarField_H
#define uniformPressureFunctionFvPatchScalarField_H

#include "dynamicPressureFvPatchScalarField.H"

#if defined( WIN32 ) || defined( WIN64 )
#include "Function1.T.H"
#else
#include "Function1.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class uniformPressureFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class uniformPressureFunctionFvPatchScalarField
:
    public dynamicPressureFvPatchScalarField
{

public:

        //- Fan flow direction
        enum fanFlowDirection
        {
            ffdIn = 0,
            ffdOut = 1
        };

        //- Fan flow directions names
        static const NamedEnum<fanFlowDirection, 2> fanFlowDirectionNames_;


private:

    // Private Data

        //- Fan curve, in function of flow rate
        autoPtr<Function1<scalar> > function_;

        //- Direction of flow through the fan relative to patch
        fanFlowDirection direction_;

        //- Name of the flux field
        const word phiName_;

        //- Volumetric flow rate minimum
        const scalar minFlowRate_;

        //- Volumetric flow rate maximum
        const scalar maxFlowRate_;

        //- Relaxation factor, because we can't have massive pressure jumps
        scalar relaxation_;

        //- Store the previous value
        scalar prevValue_;

        //- Time index of the last update
        label previousTimeIndex_;


public:

    //- Runtime type information
    TypeName("uniformPressureFunction");


    // Constructors

        //- Construct from patch and internal field
        uniformPressureFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        uniformPressureFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  uniformPressureFunctionFvPatchScalarField
        //  onto a new patch
        uniformPressureFunctionFvPatchScalarField
        (
            const uniformPressureFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new uniformPressureFunctionFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        uniformPressureFunctionFvPatchScalarField
        (
            const uniformPressureFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new uniformPressureFunctionFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
