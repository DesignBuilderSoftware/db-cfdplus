/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
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
    Foam::flowRateInletDirectionVelocityFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    Velocity inlet boundary condition either correcting the extrapolated
    velocity or creating a uniform velocity field whose direction is assumed
    to be normal to the patch, while also adjusted to match the specified flow
    rate.

    WARNING: The extrapolation option was modified and has not been tested yet.

    For a mass-based flux:
    - the flow rate should be provided in kg/s
    - if \c rho is "none" the flow rate is in m^3/s
    - otherwise \c rho should correspond to the name of the density field
    - if the density field cannot be found in the database, the user must
      specify the inlet density using the \c rhoInlet entry

    For a volumetric-based flux:
    - the flow rate is in m^3/s

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        n            | direction vector        | yes         |
        massFlowRate | mass flow rate [kg/s]   | no          |
        volumetricFlowRate | volumetric flow rate [m^3/s]| no |
        rho          | density field name      | no          | rho
        rhoInlet     | inlet density           | no          |
        extrapolateProfile | Extrapolate velocity profile | no | false
    \endtable

    Example of the boundary condition specification for a volumetric flow rate:
    \verbatim
    myPatch
    {
        type        flowRateInletDirectionVelocity;
        volumetricFlowRate  0.2;
        value       uniform (0 0 0); // placeholder
    }
    \endverbatim

    Example of the boundary condition specification for a mass flow rate with
    extrapolation:
     \verbatim
    <patchName>
    {
        type                flowRateInletDirectionVelocity;
        n                   (0.0 1.0 0.0);
        massFlowRate        0.2;
        extrapolateProfile  yes;
        rho                 rho;
        rhoInlet            1.0;
        value               uniform (0 0 0);
    }
    \endverbatim

    The \c flowRate entry is a \c Function1 of time, see Foam::Function1s.

Note
    - \c rhoInlet is required for the case of a mass flow rate, where the
      density field is not available at start-up
    - The value is positive into the domain (as an inlet)
    - May not work correctly for transonic inlets
    - Strange behaviour with potentialFoam since the U equation is not solved

See also
    Foam::fixedValueFvPatchField
    Foam::Function1s
    Foam::flowRateOutletVelocityFvPatchVectorField

SourceFiles
    flowRateInletDirectionVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef flowRateInletDirectionVelocityFvPatchVectorField_H
#define flowRateInletDirectionVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#if defined( WIN32 ) || defined( WIN64 )
#include "Function1.T.H"
#else
#include "Function1.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
           Class flowRateInletDirectionVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class flowRateInletDirectionVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Flow direction
        vectorField n_;

        //- Scaled flow direction (should be updated along with mesh updates)
        vectorField nq_;

        //- Inlet integral flow rate
        autoPtr<Function1<scalar>> flowRate_;

        //- Is volumetric?
        bool volumetric_;

        //- Name of the density field used to normalize the mass flux
        word rhoName_;

        //- Rho initialisation value (for start; if value not supplied)
        scalar rhoInlet_;

        //- Set true to extrapolate the velocity profile from the interior
        Switch extrapolateProfile_;


    // Private Member Functions

        //- Update the patch values given the appropriate density type and value
        template<class RhoType>
        void updateValues(const RhoType& rho);


public:

   //- Runtime type information
   TypeName("flowRateInletDirectionVelocity");


   // Constructors

        //- Construct from patch and internal field
        flowRateInletDirectionVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        flowRateInletDirectionVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  flowRateInletDirectionVelocityFvPatchVectorField
        //  onto a new patch
        flowRateInletDirectionVelocityFvPatchVectorField
        (
            const flowRateInletDirectionVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        flowRateInletDirectionVelocityFvPatchVectorField
        (
            const flowRateInletDirectionVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new flowRateInletDirectionVelocityFvPatchVectorField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        flowRateInletDirectionVelocityFvPatchVectorField
        (
            const flowRateInletDirectionVelocityFvPatchVectorField&,
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
                new flowRateInletDirectionVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchVectorField&, const labelList&);

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
