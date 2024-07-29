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
    Foam::incompressible::atmBLMixingLengthDissipationRateFvPatchScalarField

Group
    grpIcoRASBoundaryConditions grpInletBoundaryConditions

Description
    This boundary condition provides a turbulence dissipation, \f$\epsilon\f$
    (epsilon) inlet condition based on a specified mixing length.  The patch
    values are calculated using:

        \f[
            \epsilon_p = \frac{C_{\mu}^{0.75} k^{1.5}}{L}
        \f]

    where

    \vartable
        \epsilon_p | patch epsilon values
        C_{\mu} | Model coefficient, set to 0.09
        k       | turbulence kinetic energy
        L       | length scale
    \endvartable

    The mixing length is given by:
        \f[
            L = 100 (z/30)^(0.5)   if z > 30

            L = 100   if z <= 30
        \f]

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        z            | vertical direction      | yes         |
        phi          | flux field name         | no          | phi
        k            | turbulence kinetic energy field name | no | k
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            atmBLMixingLengthDissipationRate;
        z               (0 0 1);       //vertical direction;
        zGround         uniform 0;     // minimum z co-ordinate
        value           uniform 200;   // placeholder
    }
    \endverbatim

Note
    In the event of reverse flow, a zero-gradient condition is applied

    Recommendations for Loads on Buildings. Architectural Institute of Japan,
    http://www.aij.or.jp/jpn/symposium/2006/loads/loads.htm

See also
    Foam::inletOutletFvPatchField

SourceFiles
    atmBLMixingLengthDissipationRateFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleatmBLMixingLengthDissipationRate_H
#define incompressibleatmBLMixingLengthDissipationRate_H

#include "inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

/*---------------------------------------------------------------------------*\
 Class atmBLMixingLengthDissipationRateFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class atmBLMixingLengthDissipationRateFvPatchScalarField
:
    public inletOutletFvPatchScalarField
{
    // Private Data

        //- Direction of the z-coordinate
        vector z_;

        //- Minimum z co-ordinate
        scalarField zGround_;

        //- Name of the turbulent kinetic energy field
        word kName_;

public:

    //- Runtime type information
    TypeName("atmBLMixingLengthDissipationRate");


    // Constructors

        //- Construct from patch and internal field
        atmBLMixingLengthDissipationRateFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        atmBLMixingLengthDissipationRateFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  atmBLMixingLengthDissipationRateFvPatchScalarField
        //  onto a new patch
        atmBLMixingLengthDissipationRateFvPatchScalarField
        (
            const atmBLMixingLengthDissipationRateFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new atmBLMixingLengthDissipationRateFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Copy constructor setting internal field reference
        atmBLMixingLengthDissipationRateFvPatchScalarField
        (
            const atmBLMixingLengthDissipationRateFvPatchScalarField&,
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
                new atmBLMixingLengthDissipationRateFvPatchScalarField
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

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
