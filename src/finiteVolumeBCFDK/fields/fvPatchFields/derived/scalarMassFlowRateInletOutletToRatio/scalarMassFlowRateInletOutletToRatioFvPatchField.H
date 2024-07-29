/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's fixedProfile.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2019 FSD blueCAPE Lda  http://www.bluecape.com.pt
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
    Foam::scalarMassFlowRateInletOutletToRatioFvPatchField

Description
    This boundary condition provides a fixed value profile condition that can
    then be:
     - massFlowRate [kg/s] - the profile value is divided by the mass flow rate
       of fluid that is passing through this boundary, so that it's converted
       to mass fraction;
     - fixedMassFraction [kg/kg] - or used as-isin both mass fraction and mass
       ratio.

    Furthermore, this boundary condition needs to know which type of mass
    fraction or ratio calculation it should use:
     - massFraction - mass fraction is calculated as the mass of
       the scalar by the total mass (scalar + fluid) within the cell or face.
     - massRatio - mass ratio is calculated as the direct division between mass
       of the scalar by the mass of the fluid, within the cell or face.

Usage
    \table
        Property     | Description       | Required | Default value
        profile      | Profile function  | yes      |
        direction    | Direction of the line along which the profile is \\
                               evaluated | yes      |
        origin       | Origin of the line along which the profile is \\
                               evaluated | yes      |
        profileType  | types: massFlowRate, fixedMassFraction | yes |
        massFractionType | massFraction or massRatio | yes |
        phi          | flux field name          | no  | phi
        rho          | density field name, rhoInf or none if phi already has it
                                                | no  | rhoInf
        inletOnly    | boolean to hard-set to inlet-only | no  | false
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type       scalarMassFlowRateInletOutletToRatio;
        profile    tableFile;

        profileCoeffs
        {
            file                "UProfile";
            format              csv;
            nHeaderLine         0;
            refColumn           0;
            componentColumns    (1 2 3);
            separator           ",";
            mergeSeparators     no;
            outOfBounds         clamp;
            interpolationScheme linear;
        }
        direction        (0 1 0);
        origin           0;
        profileType      massFlowRate;
        massFractionType massFraction;
    }
    \endverbatim

    Example setting a parabolic inlet profile for the pitzDaily case:
    \verbatim
    inlet
    {
        type            scalarMassFlowRateInletOutletToRatio;

        profile         polynomial
        (
            ((1 0 0)        (0 0 0))
            ((-6200 0 0)    (2 0 0))
        );
        direction       (0 1 0);
        origin          0.0127;
        profileType     fixedMassFraction;
        massFractionType massRatio;
    }
    \endverbatim

Notes
    - The 'rhoInf' refers to using the density value defined on the
      transportProperties.

    - The profile entry is a Function1 type.  The example above gives the
      usage for supplying csv file.

See also
    Foam::fixedProfileFvPatchField
    Foam::inletOutletFvPatchField
    Foam::Function1s

SourceFiles
    scalarMassFlowRateInletOutletToRatioFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef scalarMassFlowRateInletOutletToRatioFvPatchField_H
#define scalarMassFlowRateInletOutletToRatioFvPatchField_H

#include "inletOutletFvPatchFields.H"
#include "Switch.H"

#if defined( WIN32 ) || defined( WIN64 )
#include "Function1.T.H"
#include "NamedEnum.T.H"
#else
#include "Function1.H"
#include "NamedEnum.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class scalarMassFlowRateInletOutletToRatioFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class scalarMassFlowRateInletOutletToRatioFvPatchField
:
    public inletOutletFvPatchScalarField
{
public:

    // Public Data

        //- Profile types
        enum profileTypes
        {
            massFlowRate,        //[kg/(s.m^2)]
            fixedMassFraction    //[kg/kg]
        };

        //- Profile type names
        static const NamedEnum<profileTypes, 2> profileTypeNames_;

        //- Mass fraction types
        enum massFractionTypes
        {
            massFraction,
            massRatio
        };

        //- Mass fraction type names
        static const NamedEnum<massFractionTypes, 2> massFractionTypeNames_;

private:

    // Private Data

        //- Profile function
        autoPtr<Function1<scalar>> profile_;

        //- Origin of the line along which the profile is evaluated
        scalar origin_;

        //- Direction of the line along which the profile is evaluated
        vector direction_;

        //- Profile type
        profileTypes profileType_;

        //- Mass fraction type
        massFractionTypes massFractionType_;

        //- Name of flux field
        word phiName_;

        //- Name of rho field
        word rhoName_;

        //- Name of inletOnly flag
        Switch inletOnly_;

public:

    //- Runtime type information
    TypeName("scalarMassFlowRateInletOutletToRatio");


    // Constructors

        //- Construct from patch and internal field
        scalarMassFlowRateInletOutletToRatioFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        scalarMassFlowRateInletOutletToRatioFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  scalarMassFlowRateInletOutletToRatioFvPatchField onto a new patch
        scalarMassFlowRateInletOutletToRatioFvPatchField
        (
            const scalarMassFlowRateInletOutletToRatioFvPatchField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new scalarMassFlowRateInletOutletToRatioFvPatchField(*this)
            );
        }

        //- Copy construct setting internal field reference
        scalarMassFlowRateInletOutletToRatioFvPatchField
        (
            const scalarMassFlowRateInletOutletToRatioFvPatchField&,
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
                new scalarMassFlowRateInletOutletToRatioFvPatchField(*this, iF)
            );
        }


    // Member Functions

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
