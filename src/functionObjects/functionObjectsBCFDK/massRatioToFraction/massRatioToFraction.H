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
    Copyright (C) 2013-2020 FSD blueCAPE Lda  http://www.bluecape.com.pt
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
    Foam::functionObjects::massRatioToFraction

Group
    grpFieldFunctionObjects

Description
    Calculates and writes the mass fraction, from a mass ratio.

See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    massRatioToFraction.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_massRatioToFraction_H
#define functionObjects_massRatioToFraction_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class massRatioToFraction Declaration
\*---------------------------------------------------------------------------*/

class massRatioToFraction
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the static Pressure field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("massRatioToFraction");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        massRatioToFraction
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~massRatioToFraction();

    // Member Functions

        //- Read the additional option
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
