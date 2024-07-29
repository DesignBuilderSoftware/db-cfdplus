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
    Foam::fv::volumetricHeatSource

Description
    Heat exchanger source model, in which the heat exchanger is defined as a
    selection of cells.

    It allows for two types of sources:
      1. Defining the heat power in [W]
      2. Defining the heat flux in [W/m^3]

Usage
    \table
        Property   | Description
        q          | Source value of type Function1, as a function of time
        heatSource | types are listed below
    \endtable

    Valid settings for \c heatSource:
      - power - heat power [W]
      - rate - volumetric heat [W/m^3]

    Example usage:
    \verbatim
    volumetricHeatSource1
    {
        type            volumetricHeatSource;
        active          true;

        volumetricHeatSourceCoeffs
        {
            selectionMode   cellZone;
            cellZone        heatZone;

            heatSource      rate;

            // if constant
            q               constant 10;

            // if table
            //q               table
            //(
            //  (0 0.0)
            //  (300 4.0e6)
            //);
        }
    }
    \endverbatim

SourceFiles
    volumetricHeatSource.C

\*---------------------------------------------------------------------------*/

#ifndef volumetricHeatSource_H
#define volumetricHeatSource_H

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
                    Class volumetricHeatSource Declaration
\*---------------------------------------------------------------------------*/

class volumetricHeatSource
:
    public cellSetOption
{

public:

  // Data types

        //- Enumeration listing the possible hest source input modes
        enum heatSourceType
        {
            hsPower,
            hsFlux
        };

        //- Heat source type names
        static const NamedEnum<heatSourceType, 2> heatSourceTypeNames_;

protected:

    // Protected data

        //- Heat source type
        heatSourceType heatSource_;

        //- Heat power [W] or heat volumetric flux [W/m^3]
        autoPtr<Function1<scalar>> q_;

public:

    //- Runtime type information
    TypeName("volumetricHeatSource");


    // Constructors

        //- Construct from components
        volumetricHeatSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        volumetricHeatSource(const volumetricHeatSource&) = delete;

    //- Destructor
    virtual ~volumetricHeatSource()
    {}


    // Member Functions

        // Explicit and implicit source

            //- Scalar
            virtual void addSup(fvMatrix<scalar>& eqn, const label fieldi)
            {
                NotImplemented;
            }


        // Explicit and implicit source

            //- Scalar
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            );

        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const volumetricHeatSource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
