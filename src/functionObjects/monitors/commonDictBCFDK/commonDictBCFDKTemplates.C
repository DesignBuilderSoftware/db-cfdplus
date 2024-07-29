/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Derived from OpenFOAM's surfaceFieldValue function object.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2013-2018 FSD blueCAPE Lda  http://www.bluecape.com.pt
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

\*---------------------------------------------------------------------------*/

#include "commonDictBCFDK.H"
#include "volFields.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class ReturnType, class Type>
Foam::tmp<ReturnType>
Foam::commonDictBCFDK::convertField
(
  const Type& field,
  const modeType& mode
) const
{
    switch (mode)
    {
        case mdCmptX:
        case mdCmptY:
        case mdCmptZ:
        {
            int comp = 0;
            switch (mode)
            {
                case mdCmptX:
                    comp = 0;
                    break;
                case mdCmptY:
                    comp = 1;
                    break;
                case mdCmptZ:
                    comp = 2;
                    break;
                default :
                    comp = 0;
                    break;
              }

            return field.component(comp);
            break;
        }
        case mdScalar:
        {
            return field.component(0);
            break;
        }
        case mdMag:
        default:
        {
            return mag(field);
            break;
        }
    }
}

// ************************************************************************* //
