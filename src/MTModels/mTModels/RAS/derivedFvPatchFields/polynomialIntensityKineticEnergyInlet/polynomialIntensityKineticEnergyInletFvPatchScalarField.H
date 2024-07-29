/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's atmBoundaryLayer* incompressible boundary condition
    classes.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2014-2018 FSD blueCAPE Lda  http://www.bluecape.com.pt
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
    Foam::polynomialIntensityKineticEnergyInletFvPatchScalarField

Group
    grpInletBoundaryConditions

Description
    This boundary condition provides a turbulent kinetic energy condition,
    based on user-supplied turbulence intensity, defined as a fraction of the
    mean velocity:

        \f[
            k_p = 1.5 (I |U|)^2
        \f]

    where

    \vartable
        k_p     | kinetic energy at the patch
        I       | turbulence intensity
        U       | velocity field
    \endvartable

    In the event of reverse flow, a zero-gradient condition is applied.

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        intensity    | fraction of mean field [0-1] | yes    |
        U            | velocity field name     | no          | U
        phi          | flux field name         | no          | phi
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            polynomialIntensityKineticEnergyInlet;
        z               (0 0 1); // z direction
        intensity       polynomial // y = 0.1 + 0.2x^2 + 0.0001x^3
        (
                        (0.1      0)
                        (0.2      2.0)
                        (0.0001   3.0)
        );
        zGround         0.0; //minimum z co-ordinate [m]
        value           uniform 1; // placeholder
    }
    \endverbatim

See also
    Foam::inletOutletFvPatchField

SourceFiles
    polynomialIntensityKineticEnergyInletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef polynomialIntensityKineticEnergyInletFvPatchScalarField_H
#define polynomialIntensityKineticEnergyInletFvPatchScalarField_H

#include "inletOutletFvPatchFields.H"
#if defined( WIN32 ) || defined( WIN64 )
#include "Function1.T.H"
#else
#include "Function1.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class polynomialIntensityKineticEnergyInletFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class polynomialIntensityKineticEnergyInletFvPatchScalarField
:
    public inletOutletFvPatchScalarField
{
    // Private Data

        //- Direction of the z-coordinate
        vector z_;

        //- Turbulent intensity as fraction of mean velocity
        autoPtr<Function1<scalar>> intensity_;

        //- Name of the velocity field
        word UName_;

        //- Minimum z co-ordinate
        scalarField zGround_;

public:

    //- Runtime type information
    TypeName("polynomialIntensityKineticEnergyInlet");


    // Constructors

        //- Construct from patch and internal field
        polynomialIntensityKineticEnergyInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        polynomialIntensityKineticEnergyInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  polynomialIntensityKineticEnergyInletFvPatchScalarField
        //  onto a new patch
        polynomialIntensityKineticEnergyInletFvPatchScalarField
        (
            const polynomialIntensityKineticEnergyInletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new polynomialIntensityKineticEnergyInletFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Copy constructor setting internal field reference
        polynomialIntensityKineticEnergyInletFvPatchScalarField
        (
            const polynomialIntensityKineticEnergyInletFvPatchScalarField&,
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
                new polynomialIntensityKineticEnergyInletFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchScalarField&, const labelList&);

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
