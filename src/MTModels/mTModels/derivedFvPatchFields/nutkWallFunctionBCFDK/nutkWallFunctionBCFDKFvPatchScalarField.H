/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's nutkWallFunctionFvPatchScalarField class.
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
    Foam::nutkWallFunctionBCFDKFvPatchScalarField

Description
    Turbulent viscosity wall-function boundary condition
    for high Reynolds number flows based on near-wall turbulence kinetic energy.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            nutkWallFunctionBCFDK;
        value           uniform 0;
    }
    \endverbatim

See also
    Foam::nutWallFunctionFvPatchScalarField

SourceFiles
    nutkWallFunctionBCFDKFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutkWallFunctionBCFDKFvPatchScalarField_H
#define nutkWallFunctionBCFDKFvPatchScalarField_H

#include "nutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class nutkWallFunctionBCFDKFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutkWallFunctionBCFDKFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
protected:

    // Protected Member Functions

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> nut() const;


public:

    //- Runtime type information
    TypeName("nutkWallFunctionBCFDK");


    // Constructors

        //- Construct from patch and internal field
        nutkWallFunctionBCFDKFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutkWallFunctionBCFDKFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutkWallFunctionBCFDKFvPatchScalarField
        //  onto a new patch
        nutkWallFunctionBCFDKFvPatchScalarField
        (
            const nutkWallFunctionBCFDKFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        nutkWallFunctionBCFDKFvPatchScalarField
        (
            const nutkWallFunctionBCFDKFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nutkWallFunctionBCFDKFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        nutkWallFunctionBCFDKFvPatchScalarField
        (
            const nutkWallFunctionBCFDKFvPatchScalarField&,
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
                new nutkWallFunctionBCFDKFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        //- Calculate and return the yPlus at the boundary
        virtual tmp<scalarField> yPlus() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
