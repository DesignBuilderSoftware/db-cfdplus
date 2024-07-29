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
    Foam::atmBLIntensityKineticEnergyFvPatchScalarField

Group
    grpInletBoundaryConditions

Description
    This boundary condition provides a turbulent kinetic energy condition,
    based on user-supplied inlet exponent and z direction:

        \f[
            k_p = 1.5 (I |U|)^2
        \f]

    where

    \vartable
        k_p     | kinetic energy at the patch
        I       | turbulence intensity
        U       | velocity field
    \endvartable

    The intensity is given by:
        \f[
            I = 0.1 (z/z_G)^(-a-0.05)   if z_b < z <= z_G

            I = 0.1 (z_b/z_G)^(-a-0.05)   if z <= z_b
        \f]

    In the event of reverse flow, a zero-gradient condition is applied.

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        a            | inlet exponent          | yes         |
        z            | vertical direction      | yes         |
        U            | velocity field name     | no          | U
        phi          | flux field name         | no          | phi
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type        atmBLIntensityKineticEnergy;
        a           0.22;           // inlet exponent
        z           (0 0 1);        // z direction
        zGround     uniform 0;      // minimum z co-ordinate
        value       uniform 1;      // placeholder
    }
    \endverbatim

    Recommendations for Loads on Buildings. Architectural Institute of Japan,
    http://www.aij.or.jp/jpn/symposium/2006/loads/loads.htm

See also
    Foam::inletOutletFvPatchField

SourceFiles
    atmBLIntensityKineticEnergyFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef atmBLIntensityKineticEnergyFvPatchScalarField_H
#define atmBLIntensityKineticEnergyFvPatchScalarField_H

#include "fvPatchFields.H"
#include "inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class atmBLIntensityKineticEnergyFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class atmBLIntensityKineticEnergyFvPatchScalarField
:
    public inletOutletFvPatchScalarField
{
    // Private Data

        //- Direction of the z-coordinate
        vector z_;

        //- Inlet exponent
        const scalar a_;

        //- Minimum z co-ordinate
        scalarField zGround_;

         //- Name of the velocity field
        word UName_;

public:

    //- Runtime type information
    TypeName("atmBLIntensityKineticEnergy");


    // Constructors

        //- Construct from patch and internal field
        atmBLIntensityKineticEnergyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        atmBLIntensityKineticEnergyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  atmBLIntensityKineticEnergyFvPatchScalarField
        //  onto a new patch
        atmBLIntensityKineticEnergyFvPatchScalarField
        (
            const atmBLIntensityKineticEnergyFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new atmBLIntensityKineticEnergyFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Copy constructor setting internal field reference
        atmBLIntensityKineticEnergyFvPatchScalarField
        (
            const atmBLIntensityKineticEnergyFvPatchScalarField&,
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
                new atmBLIntensityKineticEnergyFvPatchScalarField
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
