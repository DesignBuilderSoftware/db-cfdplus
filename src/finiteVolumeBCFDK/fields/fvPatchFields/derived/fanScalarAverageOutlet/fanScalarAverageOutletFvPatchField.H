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
    Foam::fanScalarAverageOutletFvPatchField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a uniform fixed value condition for the
    downstream side of a fan. It is meant to be used with the turbulent
    variables.
    The average value of the upstream side weighted by the flow rate is imposed
    on the downstream side.

    \heading Patch usage

    \table
        Property        | Description          | Required | Default value
        remotePatchName | name of remote patch | yes      |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            fanScalarAverageOutletFvPatchField;
        remotePatchName rPatch;
        phi             phi;
        value           uniform 0; // place holder
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    fanScalarAverageOutletFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef fanScalarAverageOutletFvPatchField_H
#define fanScalarAverageOutletFvPatchField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class fanScalarAverageOutletFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class fanScalarAverageOutletFvPatchField
:
    public fixedValueFvPatchScalarField
{
    // Private Data

        //- Name of the remote patch to be mapped
        word remotePatchName_;

        //- Name of the flux transporting the field
        word phiName_;


public:

    //- Runtime type information
    TypeName("fanScalarAverageOutlet");


    // Constructors

        //- Construct from patch and internal field
        fanScalarAverageOutletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fanScalarAverageOutletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fanScalarAverageOutletFvPatchField
        //  onto a new patch
        fanScalarAverageOutletFvPatchField
        (
            const fanScalarAverageOutletFvPatchField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        fanScalarAverageOutletFvPatchField
        (
            const fanScalarAverageOutletFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fanScalarAverageOutletFvPatchField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        fanScalarAverageOutletFvPatchField
        (
            const fanScalarAverageOutletFvPatchField&,
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
                new fanScalarAverageOutletFvPatchField(*this, iF)
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
