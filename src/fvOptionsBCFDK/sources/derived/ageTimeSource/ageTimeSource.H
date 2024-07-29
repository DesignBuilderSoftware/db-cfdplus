/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's buoyancyForce.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2021 FSD blueCAPE Lda  http://www.bluecape.com.pt
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
    Foam::fv::ageTimeSource

Description
    Calculates and applies the time source rho*1s to the age equation
    corresponding to the specified age field.

Usage
    Example usage:
    \verbatim
    fields          (age);
    \endverbatim

SourceFiles
    ageTimeSource.C

\*---------------------------------------------------------------------------*/

#ifndef ageTimeSource_H
#define ageTimeSource_H

#include "fvOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class ageTimeSource Declaration
\*---------------------------------------------------------------------------*/

class ageTimeSource
:
    public option
{
public:

    //- Runtime type information
    TypeName("ageTimeSource");


    // Constructors

        //- Construct from explicit source name and mesh
        ageTimeSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        ageTimeSource(const ageTimeSource&) = delete;


    // Member Functions

        // Evaluate

            //- Add explicit contribution to incompressible momentum equation
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const label fieldi
            );

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            );


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ageTimeSource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
