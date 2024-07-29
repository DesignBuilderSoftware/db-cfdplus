/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's effectivenessHeatExchangerSource.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2014-2019 FSD blueCAPE Lda  http://www.bluecape.com.pt
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
    Foam::fv::volumetricTracerSource

Description
    This utility is used to provide a source term for the scalar (tracer)
    transport equation, so it has to:
      1. Provide the tracer's units / time units, if defines a source term.
      2. Forces a fixed value as a constrain, if so defined.

    The scalar (tracer) is expected to have the units of mass fraction, so
    effectively dimensionless, hence the implementations done in function of
    the defined input value.

    It allows for three types of sources:
      1. Fixating a constant mass fraction value per unit of volume.
      2. Defining a mass flow rate per unit of volume.
      3. Defining a mass fraction flow rate per unit of volume.

    Furthermore, this volumetric condition needs to know which type of mass
    fraction or ratio calculation it should use:
     - massFraction - mass fraction is calculated as the mass of
       the scalar by the total mass (scalar + fluid) within the cell.
     - massRatio - mass ratio is calculated as the direct division between mass
       of the scalar by the mass of the fluid, within the cell.
       - This means that even though "mass fraction" is defined in the
         'sourceType', it will use mass ratio instead.

Usage
    \table
        Property   | Description
        tSource    | Source value of type Function1, as a function of time
        sourceType | types are listed below
        rho        | density field name or rhoInf
        massFractionType | massFraction or massRatio
    \endtable

    Valid settings for \c sourceType:
      - constantMassFraction - in [(kg/kg)/m^3]
      - massFlowRate - in [kg/(s.m^3)]
      - massFractionFlowRate - in [(kg/kg)/(s.m^3)]

    As for \c rho:
     - The 'rhoInf' refers to using the density value defined on the
       transportProperties.

    Example usage:
    \verbatim
    volumetricTracerSource1
    {
        type            volumetricTracerSource;
        active          true;

        volumetricTracerSourceCoeffs
        {
            selectionMode   cellZone;
            cellZone        tracerZone;

            sourceType      massFlowRate;
            rho             rhoInf;
            massFractionType massRatio;

            // if constant
            tSource         constant 10;
            tracerName      CH4;

            // if table
            //tSource         table
            //(
            //  (0 0.0)
            //  (300 4.0e6)
            //);
        }
    }
    \endverbatim

SourceFiles
    volumetricTracerSource.C

\*---------------------------------------------------------------------------*/

#ifndef volumetricTracerSource_H
#define volumetricTracerSource_H

#include "cellSetOption.H"

#if defined( WIN32 ) || defined( WIN64 )
#include "NamedEnum.T.H"
#include "Function1.T.H"
#else
#include "NamedEnum.H"
#include "Function1.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
              Class effectivenessHeatExchangerSource Declaration
\*---------------------------------------------------------------------------*/

class volumetricTracerSource
:
    public cellSetOption
{

public:

    // Public data

        //- Tracer source types
        enum sourceTypes
        {
            constantMassFraction, //[(kg/kg)/m^3]
            massFlowRate,        //[kg/(s.m^3)]
            massFractionFlowRate //[(kg/kg)/(s.m^3)]
        };

        //- Ventilator flow directions names
        static const NamedEnum<sourceTypes, 3> sourceTypeNames_;

        //- Mass fraction types
        enum massFractionTypes
        {
            massFraction,
            massRatio
        };

        //- Mass fraction type names
        static const NamedEnum<massFractionTypes, 2> massFractionTypeNames_;

protected:

    // Protected data

        //- Tracer source
        autoPtr<Function1<scalar>> tSource_;

        //- Tracer field name
        word tracerName_;

        //- Source type
        sourceTypes sourceType_;

        //- Name of rho field
        word rhoName_;

        //- Mass fraction type
        massFractionTypes massFractionType_;

private:

    // Private Member Functions

        //- Initialise heat exchanger source model
        void initialise();

public:

    //- Runtime type information
    TypeName("volumetricTracerSource");


    // Constructors

        //- Construct from components
        volumetricTracerSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        volumetricTracerSource(const volumetricTracerSource&) = delete;

    //- Destructor
    virtual ~volumetricTracerSource()
    {}


    // Member Functions

        // Public Functions

            //- Scalar
            virtual void addSup(fvMatrix<scalar>& eqn, const label fieldi);


        // Explicit and implicit source

            //- Scalar
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            );


            //- Set value on field
            virtual void constrain(fvMatrix<scalar>& eqn, const label fieldi);

        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const volumetricTracerSource&) = delete;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
