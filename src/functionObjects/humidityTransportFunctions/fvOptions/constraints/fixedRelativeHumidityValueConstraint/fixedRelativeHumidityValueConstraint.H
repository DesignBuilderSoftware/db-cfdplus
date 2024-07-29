/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's FixedValueConstraint.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2017-2018 FSD blueCAPE Lda  http://www.bluecape.com.pt
-------------------------------------------------------------------------------
License
    This file is a derivative of OpenFOAM.

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
    Foam::fv::fixedRelativeHumidityValueConstraint

Description
    Constrain the mass fraction values within a specified region, in function
    of the requested relative humidity (in percentage).

    The code was adapted from Foam::fv::FixedValueConstraint.

Usage
    For example:
    \verbatim
    humidifier
    {
        type            fixedRelativeHumidityValueConstraint;
        active          yes;
        selectionMode   cellZone;
        cellZone        porosity;

        fixedRelativeHumidityValueConstraintCoeffs
        {
            relativeHumidity  80.0; //in %
        }
    }
    \endverbatim

See also
    Foam::fvOption
    Foam::fv::FixedValueConstraint

SourceFiles
    fixedRelativeHumidityValueConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef fixedRelativeHumidityValueConstraint_H
#define fixedRelativeHumidityValueConstraint_H

#include "cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
           Class fixedRelativeHumidityValueConstraint Declaration
\*---------------------------------------------------------------------------*/

class fixedRelativeHumidityValueConstraint
:
    public cellSetOption
{
    // Private member data

        //- Relative Humidity value in fraction, therefore it's converted from
        //  percentage from the dictionary
        scalar relativeHumidity;


public:

    //- Runtime type information
    TypeName("fixedRelativeHumidityValueConstraint");


    // Constructors

        //- Construct from components
        fixedRelativeHumidityValueConstraint
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        //- Read source dictionary
        virtual bool read(const dictionary& dict);

        //- Set value on field
        virtual void constrain(fvMatrix<scalar>& eqn, const label fieldi);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
