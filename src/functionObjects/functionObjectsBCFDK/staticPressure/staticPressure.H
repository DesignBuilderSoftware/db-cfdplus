/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Derived from OpenFOAM's MachNo function object.
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

Class
    Foam::functionObjects::staticPressure

Group
    grpFieldFunctionObjects

Description
    Calculates and writes the Pressure field (rho*p) as a volScalarField.

See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    staticPressure.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_staticPressure_H
#define functionObjects_staticPressure_H

#include "fieldExpression.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class staticPressure Declaration
\*---------------------------------------------------------------------------*/

class staticPressure
:
    public fieldExpression
{
public:

    enum modeEnum
    {
        mdIsothermal,
        mdBuoyant
    };

private:

    // Private Data

        dimensionedScalar rho_;

        //- Mode names
        static const NamedEnum<modeEnum, 2> modeNames_;

        //- mode (isothermal or buoyant)
        modeEnum mode_;

    // Private Member Functions

        //- Calculate the static Pressure field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("staticPressure");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        staticPressure
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~staticPressure();


    // Member Functions

        //- Read the user-proviced data
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
