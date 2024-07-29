/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's P1 model.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2016-2018 FSD blueCAPE Lda  http://www.bluecape.com.pt
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
    Foam::radiation::P1manipulatedMedium

Description
    Works well for combustion applications where optical thickness, tau is
    large, i.e. tau = a*L > 3 (L = distance between objects)

    Assumes
     - all surfaces are diffuse
     - tends to over predict radiative fluxes from sources/sinks
       *** SOURCES NOT CURRENTLY INCLUDED ***

    The utility is the same as the OpenFOAM's except that the absorptivity and
    emissivity are no longer constant in the domain, those values can be changed
    at the boundaries. The absorptivity and emissivity are now treated as 
    boundary conditions and as so they need to be placed in 0/ folder.

SourceFiles
    P1manipulatedMedium.C

\*---------------------------------------------------------------------------*/

#ifndef radiationModelP1manipulatedMedium_H
#define radiationModelP1manipulatedMedium_H

#include "radiationModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                           Class P1manipulatedMedium Declaration
\*---------------------------------------------------------------------------*/

class P1manipulatedMedium
:
    public radiationModel
{
    // Private Data

        //- Incident radiation / [W/m^2]
        volScalarField G_;

        //- Total radiative heat flux [W/m^2]
        volScalarField qr_;

        //- Absorption coefficient
        volScalarField a_;

        //- Emission coefficient
        volScalarField e_;

        //- Emission contribution
        volScalarField E_;

public:

    //- Runtime type information
    TypeName("P1manipulatedMedium");


    // Constructors

        //- Construct from components
        P1manipulatedMedium(const volScalarField& T);

        //- Construct from components
        P1manipulatedMedium(const dictionary& dict, const volScalarField& T);

        //- Disallow default bitwise copy construction
        P1manipulatedMedium(const P1manipulatedMedium&) = delete;

    //- Destructor
    virtual ~P1manipulatedMedium();


    // Member Functions

        // Edit

            //- Solve radiation equation(s)
            void calculate();

            //- Read radiation properties dictionary
            bool read();


        // Access

            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp() const;

            //- Source term component (constant)
            virtual tmp<volScalarField::Internal> Ru() const;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const P1manipulatedMedium&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
