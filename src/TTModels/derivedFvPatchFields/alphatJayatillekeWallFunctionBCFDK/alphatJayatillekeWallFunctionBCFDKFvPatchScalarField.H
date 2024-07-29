/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's alphatJayatillekeWallFunctionFvPatchScalarField class.
    Stripped low-Reynolds handling and assume high-Reynolds only.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2021 FSD blueCAPE Lda  http://www.bluecape.com.pt
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
    Foam::compressible::alphatJayatillekeWallFunctionBCFDKFvPatchScalarField

Description
    This boundary condition provides a thermal wall function for turbulent
    thermal diffusivity (usually\c alphat) based on the Jayatilleke model.

Usage
    \table
        Property     | Description             | Required    | Default value
        Prt          | turbulent Prandtl number | no         | 0.85
        Cmu          | model coefficient       | no          | 0.09
        kappa        | Von Karman constant     | no          | 0.41
        E            | model coefficient       | no          | 9.8
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            alphatJayatillekeWallFunctionBCFDK;
        Prt             0.85;
        kappa           0.41;
        value           uniform 0; // optional value entry
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField

SourceFiles
    alphatJayatillekeWallFunctionBCFDKFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef alphatJayatillekeWallFunctionBCFDKFvPatchScalarField_H
#define alphatJayatillekeWallFunctionBCFDKFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "nutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
      Class alphatJayatillekeWallFunctionBCFDKFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class alphatJayatillekeWallFunctionBCFDKFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data

        //- Turbulent Prandtl number
        scalar Prt_;


        // Solution parameters

            static scalar maxExp_;
            static scalar tolerance_;
            static label maxIters_;


    // Private Member Functions

        //- `P' function
        scalar Psmooth(const scalar Prat) const;

        //- Calculate y+ at the edge of the thermal laminar sublayer
        scalar yPlusTherm
        (
            const nutWallFunctionFvPatchScalarField& nutw,
            const scalar P,
            const scalar Prat
        ) const;


public:

    //- Runtime type information
    TypeName("compressible::alphatJayatillekeWallFunctionBCFDK");


    // Constructors

        //- Construct from patch and internal field
        alphatJayatillekeWallFunctionBCFDKFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        alphatJayatillekeWallFunctionBCFDKFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given an
        //  alphatJayatillekeWallFunctionBCFDKFvPatchScalarField
        //  onto a new patch
        alphatJayatillekeWallFunctionBCFDKFvPatchScalarField
        (
            const alphatJayatillekeWallFunctionBCFDKFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        alphatJayatillekeWallFunctionBCFDKFvPatchScalarField
        (
            const alphatJayatillekeWallFunctionBCFDKFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new alphatJayatillekeWallFunctionBCFDKFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        alphatJayatillekeWallFunctionBCFDKFvPatchScalarField
        (
            const alphatJayatillekeWallFunctionBCFDKFvPatchScalarField&,
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
                new alphatJayatillekeWallFunctionBCFDKFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
