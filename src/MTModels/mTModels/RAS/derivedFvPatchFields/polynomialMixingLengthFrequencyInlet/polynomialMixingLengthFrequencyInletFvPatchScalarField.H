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
    Foam::incompressible::
        polynomialMixingLengthFrequencyInletFvPatchScalarField

Group
    grpIcoRASBoundaryConditions grpInletBoundaryConditions

Description
    This boundary condition provides a turbulence specific dissipation,
    \f$\omega\f$ (omega) inlet condition based on a specified mixing length.
    The patch values are calculated using:

        \f[
            \omega_p = \frac{k^{0.5}}{C_{\mu}^{0.25} L}
        \f]

    where

    \vartable
        \omega_p | patch omega values
        C_{\mu} | Model coefficient, set to 0.09
        k       | turbulence kinetic energy
        L       | length scale
    \endvartable

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        mixingLength | Length scale [m]        | yes         |
        phi          | flux field name         | no          | phi
        k            | turbulence kinetic energy field name | no | k
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            polynomialMixingLengthFrequencyInlet;
        z               (0 0 1); // z direction
        mixingLength    polynomial // y = 0.1 + 0.2x^2 + 0.0001x^3
        (
                        (0.1      0)
                        (0.2      2.0)
                        (0.0001   3.0)
        );
        zGround         0.0; //minimum z co-ordinate [m]
        value           uniform 200; // placeholder
    }
    \endverbatim

Note
    In the event of reverse flow, a zero-gradient condition is applied

See also
    Foam::inletOutletFvPatchField


SourceFiles
    polynomialMixingLengthFrequencyInletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef incompressiblepolynomialMixingLengthFrequencyInletFvPatchScalarField_H
#define incompressiblepolynomialMixingLengthFrequencyInletFvPatchScalarField_H

#include "inletOutletFvPatchFields.H"
#if defined( WIN32 ) || defined( WIN64 )
#include "Function1.T.H"
#else
#include "Function1.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

/*---------------------------------------------------------------------------*\
    Class polynomialMixingLengthFrequencyInletFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class polynomialMixingLengthFrequencyInletFvPatchScalarField
:
    public inletOutletFvPatchScalarField
{
    // Private Data

        //- Direction of the z-coordinate
        vector z_;

        //- turbulent length scale
        autoPtr<Function1<scalar> > mixingLength_;

        //- Name of the turbulent kinetic energy field
        word kName_;

        //- Minimum z co-ordinate
        scalarField zGround_;

public:

    //- Runtime type information
    TypeName("polynomialMixingLengthFrequencyInlet");


    // Constructors

        //- Construct from patch and internal field
        polynomialMixingLengthFrequencyInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        polynomialMixingLengthFrequencyInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  polynomialMixingLengthFrequencyInletFvPatchScalarField
        //  onto a new patch
        polynomialMixingLengthFrequencyInletFvPatchScalarField
        (
            const polynomialMixingLengthFrequencyInletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new polynomialMixingLengthFrequencyInletFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Copy constructor setting internal field reference
        polynomialMixingLengthFrequencyInletFvPatchScalarField
        (
            const polynomialMixingLengthFrequencyInletFvPatchScalarField&,
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
                new polynomialMixingLengthFrequencyInletFvPatchScalarField
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

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
