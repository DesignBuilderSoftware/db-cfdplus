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
    Foam::incompressible::atmBLVelocityFvPatchVectorField

Group
    grpIcoRASBoundaryConditions grpInletBoundaryConditions

Description
    This boundary condition specifies a velocity inlet profile appropriate
    for atmospheric boundary layers (ABL). Velocity is specified according to
    recommendations in the ASHRAE Handbook (chapter “Airflow around buildings”)
    co-ordinate \c z.

        \f[
            U = U_{met} (\frac{delta_met}{h_met})^a_met (\frac{h}{delta})^a
        \f]

    where
    \vartable
        U_met     | reference velocity at metereological station [m/s]
        h_met     | reference height at metereological station [m]
        delta_met | layer thickness at metereological station [m]
        a_met     | exponent at metereological station
        h         | height above ground level [m]
        delta     | layer thickness at inlet [m]
        a     | exponent at inlet

    \endvartable


    \heading Patch usage

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            atmBLVelocity;
        n               (1 0 0); //flow direction [m]
        z               (0 0 1); //vertical direction [m]
        U_met           10.0;
        h_met           10.0;
        delta_met       270;
        a_met           0.14;
        delta           370;
        a               0.22;
        zGround         0.0; //minimum z co-ordinate [m]
        value           uniform (0 0 0); // placeholder
    }
    \endverbatim

    2013 ASHRAE Handbook – Fundamentals, 2013
    (chapter “Airflow around buildings”)

SourceFiles
    atmBLVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef atmBLVelocityFvPatchVectorField_H
#define atmBLVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

/*---------------------------------------------------------------------------*\
       Class atmBLVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class atmBLVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Flow direction
        vector n_;

        //- Direction of the z-coordinate
        vector z_;

        //- Reference velocity
        const scalar Uref_;

        //- Reference height
        const scalar Href_;

        //- Reference layer thickness
        const scalar deltaref_;

        //- Reference exponent
        const scalar aref_;

        //- Inlet thickness layer
        const scalar delta_;

        //- Inlet exponent
        const scalar a_;

        //- Minimum z co-ordinate
        scalarField zGround_;

public:

    //- Runtime type information
    TypeName("atmBLVelocity");


    // Constructors

        //- Construct from patch and internal field
        atmBLVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        atmBLVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // atmBLVelocityFvPatchVectorField onto a new patch
        atmBLVelocityFvPatchVectorField
        (
            const atmBLVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new atmBLVelocityFvPatchVectorField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        atmBLVelocityFvPatchVectorField
        (
            const atmBLVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new atmBLVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Return flow direction
            const vector& n() const
            {
                return n_;
            }

            //- Return z direction
            const vector& z() const
            {
                return z_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchVectorField&, const labelList&);


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
