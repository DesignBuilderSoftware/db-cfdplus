/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's buoyantKEpsilon class. Implemented the buoyancy
    production term with the Boussinesq approximation.
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
    Foam::RASModels::buoyantkEpsilonRNGBCFDK

Group
    grpIcoRASTurbulence

Description
    Inherit the RNG k-epsilon turbulence model for incompressible flows,
    with the addition of the production term for buoyancy.

    The default model coefficients correspond to the following:
    \verbatim
        buoyantkEpsilonRNGBCFDKCoeffs
        {
            // Same as RNGkEpsilon's, plus the following:

            Cb           1.00; // It can also be 0 if the user does not want the
                               // buoyancy effects in the epsilon equation.
        }
    \endverbatim

SourceFiles
    buoyantkEpsilonRNGBCFDK.C

\*---------------------------------------------------------------------------*/

#ifndef buoyantkEpsilonRNGBCFDK_H
#define buoyantkEpsilonRNGBCFDK_H

#include "buoyantkEpsilonBaseBCFDK.H"

#if defined( WIN32 ) || defined( WIN64 )
#include "RNGkEpsilon.T.H"
#else
#include "RNGkEpsilon.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                     Class buoyantkEpsilonRNGBCFDK Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class buoyantkEpsilonRNGBCFDK
:
    public buoyantkEpsilonBaseBCFDK
    <
        RNGkEpsilon,
        dimensionedScalar,
        BasicTurbulenceModel
    >
{
protected:

    // Protected Member Functions

        //- Provide the Cmu on a type basis
        virtual dimensionedScalar Cmu() const
        {
            return this->Cmu_;
        }

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("buoyantkEpsilonRNGBCFDK");

    // Constructors

        //- Construct from components
        buoyantkEpsilonRNGBCFDK
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& type = typeName
        );

        // Disallow default bitwise copy construction
        buoyantkEpsilonRNGBCFDK(const buoyantkEpsilonRNGBCFDK&) = delete;

    //- Destructor
    virtual ~buoyantkEpsilonRNGBCFDK()
    {}

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const buoyantkEpsilonRNGBCFDK&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "buoyantkEpsilonRNGBCFDK.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
