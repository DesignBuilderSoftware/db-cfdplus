/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2015-2018 FSD blueCAPE Lda  http://www.bluecape.com.pt
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
    Foam::fanSwirlOutletVelocityFvPatchField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition will get the average flow rate from the remote
    patch and enforces a swirl, while maintaining the same average flow rate.

Usage

    \table
        Property        | Description          | Required | Default value
        remotePatchName | name of remote patch | yes      |
        phi             | flux field name      | no       | phi
        rho             | density field name   | no       | rho
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            fanSwirlOutletVelocity;
        remotePatchName rPatch;
        phi             phi;
        rho             rho;
        rpm             100;
        value           uniform (0 0 0); // place holder
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField

SourceFiles
    fanSwirlOutletVelocityFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef fanSwirlOutletVelocityFvPatchField_H
#define fanSwirlOutletVelocityFvPatchField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class fanSwirlOutletVelocityFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class fanSwirlOutletVelocityFvPatchField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Name of the remote patch to be mapped
        word remotePatchName_;

        //- Name of the flux transporting the field
        const word phiName_;

        //- Name of the density field used to normalize the mass flux
        const word rhoName_;

        //- RPM
        const scalar rpm_;


public:

    //- Runtime type information
    TypeName("fanSwirlOutletVelocity");


    // Constructors

        //- Construct from patch and internal field
        fanSwirlOutletVelocityFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fanSwirlOutletVelocityFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fanSwirlOutletVelocityFvPatchField
        //  onto a new patch
        fanSwirlOutletVelocityFvPatchField
        (
            const fanSwirlOutletVelocityFvPatchField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector>> clone() const
        {
            return tmp<fvPatchField<vector>>
            (
                new fanSwirlOutletVelocityFvPatchField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        fanSwirlOutletVelocityFvPatchField
        (
            const fanSwirlOutletVelocityFvPatchField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector>> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector>>
            (
                new fanSwirlOutletVelocityFvPatchField(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Name of the remote patch to be mapped
            const word& remotePatchName() const
            {
                return remotePatchName_;
            }

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
