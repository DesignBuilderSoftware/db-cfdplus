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
    Foam::functionObjects::MachNoBCFDK

Group
    grpFieldFunctionObjects

Description
    Calculates and writes the Mach number as a volScalarField.

See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    MachNoBCFDK.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_MachNoBCFDK_H
#define functionObjects_MachNoBCFDK_H

#include "fieldExpression.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class MachNoBCFDK Declaration
\*---------------------------------------------------------------------------*/

class MachNoBCFDK
:
    public fieldExpression
{
    // Private Data

        //- Switch on whether it's isothermal or not
        Switch isIsothermal_;

        dimensionedScalar TRef_;
        dimensionedScalar gamma_;
        dimensionedScalar gasConst_;
        dimensionedScalar airMass_;

    // Private Member Functions

        //- Calculate the Mach number field in function of temperature and
        //  return true if successful
        bool compressibleCalc();
        template<class TempType> bool incompressibleCalc(const TempType &temp);

        //- Calculate the Mach number field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("MachNoBCFDK");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        MachNoBCFDK
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~MachNoBCFDK();


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
