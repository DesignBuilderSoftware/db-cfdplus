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
    Foam::fanRedirectedOutletVelocityFvPatchField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition will get the average flow rate from the remote
    patch and enforces a new direction, while maintaining the same average flow
    rate.

    \heading Patch usage

    \table
        Property        | Description          | Required | Default value
        remotePatchName | name of remote patch | yes      |
        U               | velocity field name  | no       | U
        phi             | flux field name      | no       | none
        rho             | density field name   | no       | none
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            fanRedirectedOutletVelocity;
        remotePatchName rPatch;
        U               U;
        phi             none;
        rho             none;
        direction       uniform (1 0 0);
        value           uniform (0 0 0); // place holder
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    fanRedirectedOutletVelocityFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef fanRedirectedOutletVelocityFvPatchField_H
#define fanRedirectedOutletVelocityFvPatchField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class fanRedirectedOutletVelocityFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class fanRedirectedOutletVelocityFvPatchField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Name of the remote patch to be mapped
        const word remotePatchName_;

        //- Name of the velocity field
        const word UName_;

        //- Name of the flux field
        const word phiName_;

        //- Name of the density field
        const word rhoName_;

        //- Allow controlling the fixesValue method
        const Switch autoAdjustable_;

        //- Inlet direction taking into account normal of patch
        vectorField inletDir_;

        //- Inlet direction given by the user
        vectorField direction_;


public:

    //- Runtime type information
    TypeName("fanRedirectedOutletVelocity");


    // Constructors

        //- Construct from patch and internal field
        fanRedirectedOutletVelocityFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fanRedirectedOutletVelocityFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fanRedirectedOutletVelocityFvPatchField
        //  onto a new patch
        fanRedirectedOutletVelocityFvPatchField
        (
            const fanRedirectedOutletVelocityFvPatchField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        fanRedirectedOutletVelocityFvPatchField
        (
            const fanRedirectedOutletVelocityFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new fanRedirectedOutletVelocityFvPatchField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        fanRedirectedOutletVelocityFvPatchField
        (
            const fanRedirectedOutletVelocityFvPatchField&,
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
                new fanRedirectedOutletVelocityFvPatchField(*this, iF)
            );
        }


    // Member Functions

        // Attributes

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return !autoAdjustable_;
            }

        // Access

            //- Name of the remote patch to be mapped
            const word& remotePatchName() const
            {
                return remotePatchName_;
            }

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
