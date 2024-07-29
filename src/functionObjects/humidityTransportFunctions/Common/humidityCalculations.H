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

Description
    Provides template functions that are used by the solver, boundary
    conditions and fvOptions.

SourceFiles
    humidityCalculations.C

\*---------------------------------------------------------------------------*/

#ifndef humidityCalculations_H
#define humidityCalculations_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    template<class Type>
    tmp<Type> calcSaturationPressure(const Type &dimlessTemperature);

    //Note: mixturePressure must be the total pressure in Pa
    template<typename Type1, typename Type2>
    void calcHumidityFraction
    (
        Type2 &result,
        const Type1 &Rfraction,
        const Type2 &relativeHumidity,
        const Type2 &saturationPressure,
        const Type2 &mixturePressure
    );

    template<typename Type>
    void calcMassFraction(Type &result, const Type &humidityFraction);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "humidityCalculations.C"
#endif

#endif

// ************************************************************************* //
