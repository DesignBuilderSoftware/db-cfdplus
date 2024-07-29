/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::RASModels::buoyantkOmegaSSTBCFDK

Group
    grpIcoRASTurbulence

Description
    Standard k-epsilon turbulence model for incompressible flows, with the
    addition of the production term for buoyancy.

    The default model coefficients correspond to the following:
    \verbatim
        buoyantkOmegaSSTBCFDKCoeffs
        {
            // Same as kOmegaSST's, plus the following:

            Cb           1.00; // It can also be 0 if the user does not want
                               // the buoyancy effects in the omega equation.
        }
    \endverbatim

    Followed documentation in subsection 2.2.2 of the CFX theory manual, ANSYS
    16.2.3.

SourceFiles
    buoyantkOmegaSSTBCFDK.C

\*---------------------------------------------------------------------------*/

#ifndef buoyantkOmegaSSTBCFDK_H
#define buoyantkOmegaSSTBCFDK_H

#include "kOmegaSST.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                     Class buoyantkOmegaSSTBCFDK Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class buoyantkOmegaSSTBCFDK
:
    public kOmegaSST<BasicTurbulenceModel>
{

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cb_;

    // Protected Member Functions

        tmp<volScalarField> Gcoef() const;

        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> omegaSource() const;

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("buoyantkOmegaSSTBCFDK");

    // Constructors

        //- Construct from components
        buoyantkOmegaSSTBCFDK
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& type = typeName
        );

        // Disallow default bitwise copy construction and assignment
        buoyantkOmegaSSTBCFDK(const buoyantkOmegaSSTBCFDK&) = delete;

    //- Destructor
    virtual ~buoyantkOmegaSSTBCFDK()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const buoyantkOmegaSSTBCFDK&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "buoyantkOmegaSSTBCFDK.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
