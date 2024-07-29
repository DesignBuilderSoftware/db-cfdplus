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
    Foam::RASModels::buoyantkEpsilonBaseBCFDK

Group
    grpIcoRASTurbulence

Description
    Template k-epsilon turbulence model for adding the production term for
    buoyancy.

    The default model coefficients correspond to the following:
    \verbatim
        buoyantkEpsilonBaseBCFDKCoeffs
        {
            // Same as kEpsilon's, plus the following defaults:

            Cbe         1.00; // It can also be 0 if the user does not want the
                              // buoyancy effects in the epsilon equation.
            Cbk         1.00; // It can also be 0 if the user does not want the
                              // buoyancy effects in the k equation.

            Tanhe       false; // when true, uses tanh(u/v) to be a function of
                               // gravitationally aligned velocity in the
                               // epsilon source term
            Tanhk       false; // when true, uses tanh(u/v) to be a function of
                               // gravitationally aligned velocity in the
                               // epsilon source term
        }
    \endverbatim

    where \c buoyantkEpsilonBaseBCFDKCoeffs should be renamed to the derived
    class instead.

See also
    Foam::RASModels::kEpsilon

SourceFiles
    buoyantkEpsilonBaseBCFDK.C

\*---------------------------------------------------------------------------*/

#ifndef buoyantkEpsilonBaseBCFDK_H
#define buoyantkEpsilonBaseBCFDK_H

#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class buoyantkEpsilonBaseBCFDK Declaration
\*---------------------------------------------------------------------------*/

template
<
    template<class> class kEpsilonBaseModel,
    class cmuType,
    class BasicTurbulenceModel
>
class buoyantkEpsilonBaseBCFDK
:
    public kEpsilonBaseModel<BasicTurbulenceModel>
{

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cbe_;
            dimensionedScalar Cbk_;
            Switch Tanhe_;
            Switch Tanhk_;

    // Protected Member Functions

        tmp<volScalarField> tanhU(
            const uniformDimensionedVectorField& g
            ) const;

        tmp<volScalarField> Gcoef(
            const uniformDimensionedVectorField& g
            ) const;

        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> epsilonSource() const;

        //- Provide the Cmu on a type basis
        virtual cmuType Cmu() const
        {
            NotImplemented;
            return cmuType();
        }

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information - not implemented here
    // TypeName("buoyantkEpsilonBaseBCFDK");

    // Constructors

        //- Construct from components
        buoyantkEpsilonBaseBCFDK
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& type = word() // typeName
        );

        // Disallow default bitwise copy construction
        buoyantkEpsilonBaseBCFDK(const buoyantkEpsilonBaseBCFDK&) = delete;

    //- Destructor
    virtual ~buoyantkEpsilonBaseBCFDK()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const buoyantkEpsilonBaseBCFDK&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "buoyantkEpsilonBaseBCFDK.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
