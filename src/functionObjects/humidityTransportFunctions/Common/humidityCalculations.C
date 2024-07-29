/*---------------------------------------------------------------------------*\
    Copyright (C) 2017 FSD blueCAPE Lda  http://www.bluecape.com.pt
-------------------------------------------------------------------------------
License
    This file is derivative of OpenFOAM.

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
tmp<Type> calcSaturationPressure
(
    const Type &dimlessTemperature
)
{

    /*
    Source: ASHRAE Fundamentals 2009 Handbook (from
    https://openfoamwiki.net/index.php/Contrib_massBuoyantBoussinesqSimpleFoam)
    pSaturation =
      EXP
      (
        -0.58002206*POWER(10,4)*POWER(T,-1)
        +0.13914993*10
        -0.48640239*POWER(10,-1)*T
        +0.41764768*POWER(10,-4)*POWER(T,2)
        -0.14452093*POWER(10,-7)*POWER(T,3)
        +0.65459673*10*LN(T)
      )
    */

    return tmp<Type>
    (
        new Type
        (
            exp
            (
                -0.58002206*10000.0*pow(dimlessTemperature,-1)
              + 0.13914993*10.0
              - 0.48640239*0.1*dimlessTemperature
              + 0.41764768*0.0001*pow(dimlessTemperature,2)
              - 0.14452093*0.0000001*pow(dimlessTemperature,3)
              + 0.65459673*10.0*log(dimlessTemperature)
            )
        )
    );
}


template<typename Type1, typename Type2>
void calcHumidityFraction
(
    Type2 &result,
    const Type1 &Rfraction,
    const Type2 &relativeHumidity,
    const Type2 &saturationPressure,
    const Type2 &mixturePressure
)
{
    result ==
        Rfraction * relativeHumidity * saturationPressure
      / (mixturePressure - relativeHumidity*saturationPressure);
}


template<typename Type>
void calcMassFraction(Type &result, const Type &humidityFraction)
{
  result == humidityFraction / (scalar(1) + humidityFraction);
}

} // End namespace Foam

// ************************************************************************* //
