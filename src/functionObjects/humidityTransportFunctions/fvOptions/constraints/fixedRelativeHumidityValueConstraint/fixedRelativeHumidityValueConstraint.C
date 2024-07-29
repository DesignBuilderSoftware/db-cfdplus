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

\*---------------------------------------------------------------------------*/

#include "fixedRelativeHumidityValueConstraint.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "Common/humidityCalculations.H"

#if defined( WIN32 ) || defined( WIN64 )
#include "DimensionedField.T.H"
#else
#include "DimensionedField.H"
#endif

// * * * * * * * * * * * * * * * Local Classes * * * * * * * * * * * * * * * //

namespace Foam
{

    class scalarFieldOL : public scalarField
    {
      public:
        scalarFieldOL(const label size) :
            scalarField(size)
        {}

        scalarFieldOL(const label size, const scalar& t) :
            scalarField(size, t)
        {}

        scalarFieldOL(const tmp<scalarField>& tf) :
            scalarField(tf)
        {}

        void operator==(const scalarFieldOL& other)
        {
            scalarField::operator=(other);
        }
    };

};

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(fixedRelativeHumidityValueConstraint, 0);
    addToRunTimeSelectionTable
    (
        option,
        fixedRelativeHumidityValueConstraint,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::fixedRelativeHumidityValueConstraint::
fixedRelativeHumidityValueConstraint
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::fixedRelativeHumidityValueConstraint::read
(
    const dictionary& dict
)
{
    if (cellSetOption::read(dict))
    {
        relativeHumidity =
            coeffs_.lookup<scalar>("relativeHumidity") / scalar(100.0);

        fieldNames_.setSize(1, "x_H2O");
        applied_.setSize(1, false);

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::fv::fixedRelativeHumidityValueConstraint::constrain
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    DebugInfo
            << "fixedRelativeHumidityValueConstraint"
            << "::constrain for source " << name_ << endl;

    //Need to convert the values:
    // - from relative humidity fraction to humidity fraction
    // - and then from humidity fraction to mass fraction

    //Load in the list of cells with the user-defined relative humidity
    scalarFieldOL phi_H2O(cells_.size(), relativeHumidity);

    //Retrieve RAirOverRH2O from transportProperties
    const dictionary& transportProperties =
        mesh_.lookupObject<IOdictionary>("transportProperties");
    const dimensionedScalar RAirOverRH2Odimed
    (
        transportProperties.lookup("RAirOverRH2O")
    );
    const scalar RAirOverRH2O = RAirOverRH2Odimed.value();

    //Retrieve saturationPressure and airPressure and transfer the respective
    //cells to temporary list
    const volScalarField& saturationPressureField =
        mesh_.lookupObject<volScalarField>("saturationPressure");
    const volScalarField& airPressureField =
        mesh_.lookupObject<volScalarField>("airPressure");

    scalarFieldOL saturationPressure(cells_.size());
    scalarFieldOL airPressure(cells_.size());
    label locali = 0;
    forAll(cells_, celli)
    {
        saturationPressure[locali] = saturationPressureField[celli];
        airPressure[locali] = airPressureField[celli];
        locali++;
    }

    //Create the w_H2O and x_H2O fields
    scalarFieldOL w_H2O(cells_.size());
    scalarFieldOL x_H2O(cells_.size());

    //calculate w_H2O
    calcHumidityFraction
    (
        w_H2O,
        RAirOverRH2O,
        phi_H2O,
        saturationPressure,
        airPressure
    );

    //calculate x_H2O
    calcMassFraction(x_H2O, w_H2O);

    eqn.setValues(cells_, x_H2O);
}


// ************************************************************************* //
