/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Modified and developed for blueCFD(R)-Kernel
    Derived from blueCFD(R)-Kernel's massAccounting function object.
    Copyright (C) 2020 FSD blueCAPE Lda  http://www.bluecape.com.pt
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
    Foam::functionObjects::resetBuoyantPressures

Group
    grpFieldFunctionObjects

Description
    Rebalance the p pressure in function of p_rgh at write time, effectively
    reverting back p_rgh and setting p instead.
    Warning: Should only be used at time step 0, as a pre-processing step, by
    calling the solver in postProcess mode.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    resetBuoyantPressures.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_resetBuoyantPressures_H
#define functionObjects_resetBuoyantPressures_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class resetBuoyantPressures Declaration
\*---------------------------------------------------------------------------*/

class resetBuoyantPressures
:
    public fvMeshFunctionObject
{
    // Private Member Functions

        //- Rebalance the p pressure in function of p_rgh
        void rebalance();


public:

    //- Runtime type information
    TypeName("resetBuoyantPressures");


    // Constructors

        //- Construct from Time and dictionary
        resetBuoyantPressures
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        resetBuoyantPressures(const resetBuoyantPressures&) = delete;

    //- Destructor
    virtual ~resetBuoyantPressures();


    // Member Functions

        //- Read the resetBuoyantPressures data
        virtual bool read(const dictionary&);

        //- Do not do anything
        virtual bool execute();

        //- Calculate the resetBuoyantPressures
        virtual bool write();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const resetBuoyantPressures&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
