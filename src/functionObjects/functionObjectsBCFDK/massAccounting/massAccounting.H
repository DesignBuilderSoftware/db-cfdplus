/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Derived from OpenFOAM's scalarTransport function object.
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
    Foam::functionObjects::massAccounting

Group
    grpFieldFunctionObjects

Description
    Checks how the mass balance and total inlet is currently going and outputs
    it only to standard output.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    massAccounting.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_massAccounting_H
#define functionObjects_massAccounting_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class massAccounting Declaration
\*---------------------------------------------------------------------------*/

class massAccounting
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Name of flux field (optional)
        word phiName_;


    // Private Member Functions

        //- Report the mass balances
        void report();


public:

    //- Runtime type information
    TypeName("massAccounting");


    // Constructors

        //- Construct from Time and dictionary
        massAccounting
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        massAccounting(const massAccounting&) = delete;

    //- Destructor
    virtual ~massAccounting();


    // Member Functions

        //- Read the massAccounting data
        virtual bool read(const dictionary&);

        //- Do not do anything
        virtual bool execute();

        //- Calculate the massAccounting
        virtual bool write();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const massAccounting&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
