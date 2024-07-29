/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's externalWallHeatFluxTemperatureFvPatchScalarField class
    Reimplemented to calculate value instead of gradient.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2023 FSD blueCAPE Lda  http://www.bluecape.com.pt
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
    Foam::wallHeatFluxTemperatureBCFDKFvPatchScalarField

Description
    This boundary condition estimates the temperature on a boundary wall,
    in function of an estimated energy quantity being transferred, in one
    of two heat modes:

      - fixed power: supply Q
      - fixed heat flux: supply q

    where:
    \vartable
        Q  | Power [W]
        q  | Heat flux [W/m^2]
    \endvartable

    It re-calculates the actual temperature value in function of the chosen
    heat transfer coefficient mode \c htcMode:

      - 'zeroGradient': use zero gradient
      - 'defined': based on a heat transfer value \c h
      - 'fluidKappa': use the fluid's heat transfer coefficient

Usage
    \table
    Property     | Description                 | Required | Default value
    mode         | 'power' or 'flux' | yes |
    Q            | Power [W]                   | for mode 'power'     |
    q            | Heat flux [W/m^2]           | for mode 'flux'     |
    htcMode      | 'zeroGradient', 'defined' or 'fluidKappa' |
    h         | Heat transfer coefficient [W/(m^2.K)] | for htcMode 'defined' |
    Tmin       | If defined, limit the minimum temperature  |  no      |
    Tmax       | If defined, limit the maximum temperature  |  no      |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            wallHeatFluxTemperatureBCFDK;
        mode            flux;
        q               uniform 10.0;

        mode            power;
        Q               10.0;

        htcMode         defined;
        h               1.5;

        Tmin     260;
        Tmax     500;

        value           $internalField;
    }
    \endverbatim

See also
    Foam::temperatureCoupledBase
    Foam::fixedValueFvPatchScalarField

SourceFiles
    wallHeatFluxTemperatureBCFDKFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef wallHeatFluxTemperatureBCFDKFvPatchScalarField_H
#define wallHeatFluxTemperatureBCFDKFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "temperatureCoupledBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class wallHeatFluxTemperatureBCFDKFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class wallHeatFluxTemperatureBCFDKFvPatchScalarField
:
    public fixedValueFvPatchScalarField,
    public temperatureCoupledBase
{
public:

    // Public data

        //- Operation mode enumeration
        enum operationMode
        {
            fixedPower,
            fixedHeatFlux
        };

        static const NamedEnum<operationMode, 2> operationModeNames;

        //- Operation htcMode enumeration
        enum operationHtcMode
        {
            zeroGradient,
            defined,
            fluidKappa
        };

        static const NamedEnum<operationHtcMode, 3> htcModeNames;


private:

    // Private Data

        //- Operation mode
        operationMode mode_;

        //- Heat transfer coefficient mode
        operationHtcMode htcMode_;

        //- Heat power [W]
        scalar Q_;

        //- Heat flux [W/m^2]
        scalarField q_;

        //- Heat transfer coefficient [W/(m^2.K)]
        scalar h_;

        //- Temperature min-max values [K]
        scalar Tmin_;
        scalar Tmax_;

public:

    //- Runtime type information
    TypeName("wallHeatFluxTemperatureBCFDK");


    // Constructors

        //- Construct from patch and internal field
        wallHeatFluxTemperatureBCFDKFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        wallHeatFluxTemperatureBCFDKFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  wallHeatFluxTemperatureBCFDKFvPatchScalarField
        //  onto a new patch
        wallHeatFluxTemperatureBCFDKFvPatchScalarField
        (
            const wallHeatFluxTemperatureBCFDKFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        wallHeatFluxTemperatureBCFDKFvPatchScalarField
        (
            const wallHeatFluxTemperatureBCFDKFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new wallHeatFluxTemperatureBCFDKFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        wallHeatFluxTemperatureBCFDKFvPatchScalarField
        (
            const wallHeatFluxTemperatureBCFDKFvPatchScalarField&,
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
                new wallHeatFluxTemperatureBCFDKFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Do not allow manipulation of the boundary values
            virtual bool fixesValue() const
            {
                return true;
            }


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
            void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
