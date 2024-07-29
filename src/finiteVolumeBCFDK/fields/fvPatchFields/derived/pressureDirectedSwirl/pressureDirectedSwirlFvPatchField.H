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
    Foam::pressureDirectedSwirlFvPatchField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition will get the average flow rate from the current
    patch and enforces a swirl, while maintaining the same average flow rate.

    Flow direction is always needed, while either RPM or diffuser vector are
    alternatives, not both can be activated.

Usage

    \table
        Property        | Description          | Required | Default value
        phi             | flux field name      | no       | phi
        rho             | density field name   | no       | rho
        direction       | direction of flow through patch [in/out] | yes |
        rpm             | rotations per minute   | no       | 0.0
        diffuser        | diffuser vector direction | no | <not used>
        minFlowRate     | minimum flow-rate allowed (>=0) | no | not limited
        maxFlowRate     | maximum flow-rate allowed (>0) | no | not limited
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatchSwirl
    {
        type            pressureDirectedSwirl;
        phi             phi;
        rho             rho;
        rpm             100;
        direction       in;
        value           uniform (0 0 0); // place holder
    }

    myPatchNoSwirl
    {
        type            pressureDirectedSwirl;
        phi             phi;
        rho             rho;
        direction       in;
        diffuser        uniform (1 0 0);
        value           uniform (0 0 0); // place holder
    }
    \endverbatim

Note
    The function entry is a DataEntry type, able to describe time
    varying functions. The example above gives the usage for supplying a
    constant value.

    This class was based on OpenFOAM's fanPressureFvPatchField and based on
    blueCFD-Kernel's class fanSwirlOutletVelocityFvPatchField. Later with
    OpenFOAM 8, adopted pressureDirectedInletVelocity as a basis.

See also
    Foam::fixedValueFvPatchField
    Foam::fanSwirlOutletVelocityFvPatchField (FSD blueCAPE)

SourceFiles
    pressureDirectedSwirlFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef pressureDirectedSwirlFvPatchField_H
#define pressureDirectedSwirlFvPatchField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class pressureDirectedSwirlFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class pressureDirectedSwirlFvPatchField
:
    public fixedValueFvPatchVectorField
{
public:

        //- Ventilator flow direction
        enum ventilatorFlowDirection
        {
            vfdIn,
            vfdOut
        };

        //- Ventilator flow directions names
        static const NamedEnum<ventilatorFlowDirection, 2>
            ventilatorFlowDirectionNames_;


private:

    // Private Data

        //- Name of the flux transporting the field
        const word phiName_;

        //- Name of the density field used to normalize the mass flux
        const word rhoName_;

        //- RPM
        const scalar rpm_;

        //- RPM is set if true, otherwise diffuser mode is used
        const bool rpmOrDiffuser_;

        //- Volumetric flow rate minimum
        const scalar minFlowRate_;

        //- Volumetric flow rate maximum
        const scalar maxFlowRate_;

        //- Diffuser direction
        vectorField diffuserDir_;

        //- Direction of flow through the ventilator relative to patch
        ventilatorFlowDirection direction_;

public:

    //- Runtime type information
    TypeName("pressureDirectedSwirl");


    // Constructors

        //- Construct from patch and internal field
        pressureDirectedSwirlFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pressureDirectedSwirlFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given pressureDirectedSwirlFvPatchField
        //  onto a new patch
        pressureDirectedSwirlFvPatchField
        (
            const pressureDirectedSwirlFvPatchField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new pressureDirectedSwirlFvPatchField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        pressureDirectedSwirlFvPatchField
        (
            const pressureDirectedSwirlFvPatchField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new pressureDirectedSwirlFvPatchField(*this, iF)
            );
        }


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchVectorField&, const labelList&);

        // Evaluation functions

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
