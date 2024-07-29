/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Somewhat based on OpenFOAM's old turbulentHeatFluxTemperature, now named
        externalWallHeatFluxTemperature.
    Mixed boundary condition mode based on inletOutlet/outletInlet.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2014-2021 FSD blueCAPE Lda  http://www.bluecape.com.pt
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
    Foam::scalarFluxFickFirstLawFvPatchScalarFieldBCFDK

Description
    This utility is used at a walls for the tracer scalar transport equation
    (carried in mass fraction kg/kg) to assign a scalar flux, which may be only
    defined as a flux definition of mass flux rate in [kg/(s.m^2)].

    Warning: We did try to implement with mass fraction flux rate in
    [(kg/kg)/(s.m^2)], but it simply did not make sense, since the wall cannot
    inject the main fluid.

    The defined flux is then converted to be transported into/from the domain
    through the gradient condition and Fick's first law of diffusion:

      J = -D . \rho_{fluid} . (d \Phi) / (d x)

    Where for mass fraction, it becomes:
      - J is the diffusion flux, in [kg/(m^2.s)]
      - D is the diffusion coefficient or diffusivity, in [m^2/s]
      - \rho_{fluid} is the density of the base transport fluid or of the
        mixture, in [kg/m^3]
      - \Phi is the mass fraction, in [kg/kg]
      - x is position, in [m]

    The gradient is limited to try to never impose a negative value at the
    boundary, in order to avoid hindering the scalar field.

    Furthermore, this boundary condition needs to know which type of mass
    fraction or ratio calculation it should use:
     - massFraction - mass fraction is calculated as the mass of
       the scalar by the total mass (scalar + fluid) within the cell or face.
     - massRatio - mass ratio is calculated as the direct division between mass
       of the scalar by the mass of the fluid, within the cell or face.

Usage

    \table
        Property   | Description                     | Required | Default value
        molDiffEff | effective molecular diffusivity | yes      | DT_scalarName
        flux       | flux value in  [kg/(s.m^2)]      | yes      | 0
        rho        | density field name or rhoInf | no  | rhoInf
        massFractionType | massFraction or massRatio | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            scalarFluxFickFirstLawBCFDK;
        flux            uniform 1e-5;
        molDiffEff      DT_CH4;
        massFractionType massFraction;
        value           $internalField; // place holder
    }
    \endverbatim

Notes
    - The 'rhoInf' refers to using the density value defined on the
      transportProperties.

    - The units of the field indicated by \c molDiffEff are [m^2/s]

SourceFiles
    scalarFluxFickFirstLawFvPatchScalarFieldBCFDK.C

\*---------------------------------------------------------------------------*/

#ifndef scalarFluxFickFirstLawFvPatchScalarFieldBCFDK_H
#define scalarFluxFickFirstLawFvPatchScalarFieldBCFDK_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"

#if defined( WIN32 ) || defined( WIN64 )
#include "NamedEnum.T.H"
#else
#include "NamedEnum.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class scalarFluxFickFirstLawFvPatchScalarFieldBCFDK Declaration
\*---------------------------------------------------------------------------*/

class scalarFluxFickFirstLawFvPatchScalarFieldBCFDK
:
    public mixedFvPatchScalarField
{
public:

    // Public data

        //- Mass fraction types
        enum massFractionTypes
        {
            massFraction,
            massRatio
        };

        //- Mass fraction type names
        static const NamedEnum<massFractionTypes, 2> massFractionTypeNames_;

private:
    // Private Data

        //- Flux value, in either
        scalarField fluxValue_;

        //- Name of effective thermal diffusivity field
        word molDiffEffName_;

        //- Name of rho field
        word rhoName_;

        //- Mass fraction type
        massFractionTypes massFractionType_;

public:

    //- Runtime type information
    TypeName("scalarFluxFickFirstLawBCFDK");


    // Constructors

        //- Construct from patch and internal field
        scalarFluxFickFirstLawFvPatchScalarFieldBCFDK
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        scalarFluxFickFirstLawFvPatchScalarFieldBCFDK
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  scalarFluxFickFirstLawFvPatchScalarFieldBCFDK onto
        //  a new patch
        scalarFluxFickFirstLawFvPatchScalarFieldBCFDK
        (
            const scalarFluxFickFirstLawFvPatchScalarFieldBCFDK&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        scalarFluxFickFirstLawFvPatchScalarFieldBCFDK
        (
            const scalarFluxFickFirstLawFvPatchScalarFieldBCFDK&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new scalarFluxFickFirstLawFvPatchScalarFieldBCFDK
                (
                    *this
                )
            );
        }

        //- Copy constructor setting internal field reference
        scalarFluxFickFirstLawFvPatchScalarFieldBCFDK
        (
            const scalarFluxFickFirstLawFvPatchScalarFieldBCFDK&,
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
                new scalarFluxFickFirstLawFvPatchScalarFieldBCFDK
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


        // I-O

            //- Write
            virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
