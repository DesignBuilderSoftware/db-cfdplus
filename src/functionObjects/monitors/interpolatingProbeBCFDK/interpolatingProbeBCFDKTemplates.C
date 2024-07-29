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

#include "interpolatingProbeBCFDK.H"
#include "volFields.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::interpolatingProbesBCFDK::probePosition
(
    const word& fieldName,
    const modeType& mode,
    const word& identifier
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;
    const bool ok = obr_.foundObject<fieldType>(fieldName);

    if (ok)
    {
        const label proci = Pstream::myProcNo();

        const fieldType& origField = obr_.lookupObject<fieldType>(fieldName);
        const volScalarField field
        (
            convertField<volScalarField, fieldType>(origField, mode)
        );

        const scalar unsetVal(-VGREAT);

        List<scalar> sampledValues(Pstream::nProcs());
        sampledValues[proci] = unsetVal;

        // Must always start up the interpolator, which needs parallel
        // capabilities
        autoPtr<interpolation<scalar>> interpolator
        (
            interpolation<scalar>::New(interpolationScheme_, field)
        );

        if (faceProbe_ >= 0)
        {
            sampledValues[proci] = interpolator().interpolate
            (
                probePosition_,
                elementProbe_,
                faceProbe_
            );
        }

        Pstream::gatherList(sampledValues);

        const label maxI = findMax(sampledValues);

        Log
            << identifier
            << token::ASSIGN << token::SPACE
            << sampledValues[maxI] << token::END_STATEMENT
            << token::SPACE;

        if (logToFile_)
        {
            file() << tab << sampledValues[maxI];
        }
    }

    return ok;
}

// ************************************************************************* //
