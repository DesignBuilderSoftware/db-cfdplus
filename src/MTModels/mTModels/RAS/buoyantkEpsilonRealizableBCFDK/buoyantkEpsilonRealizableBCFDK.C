/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's kEpsilon model. Added the buyoancy production term.
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

\*---------------------------------------------------------------------------*/

#include "buoyantkEpsilonRealizableBCFDK.H"
#include "uniformDimensionedFields.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> buoyantkEpsilonRealizableBCFDK<BasicTurbulenceModel>
    ::Cmu() const
{
    //Code lifted from 'realizableKE::correctNut()'
    tmp<volTensorField> tgradU = fvc::grad(this->U_);
    volScalarField S2(2*magSqr(dev(symm(tgradU()))));
    volScalarField magS(sqrt(S2));
    return reimplemented_rCmu(tgradU(), S2, magS);
}


template<class BasicTurbulenceModel>
tmp<volScalarField> buoyantkEpsilonRealizableBCFDK<BasicTurbulenceModel>
    ::reimplemented_rCmu
(
    const volTensorField& gradU,
    const volScalarField& S2,
    const volScalarField& magS
) const
{
    tmp<volSymmTensorField> tS = dev(symm(gradU));
    const volSymmTensorField& S = tS();

    volScalarField W
    (
        (2*sqrt(2.0))*((S&S)&&S)
       /(
            magS*S2
          + dimensionedScalar(dimensionSet(0, 0, -3, 0, 0), small)
        )
    );

    tS.clear();

    volScalarField phis
    (
        (1.0/3.0)*acos(min(max(sqrt(6.0)*W, -scalar(1)), scalar(1)))
    );
    volScalarField As(sqrt(6.0)*cos(phis));
    volScalarField Us(sqrt(S2/2.0 + magSqr(skew(gradU))));

    return 1.0/(this->A0_ + As*Us*this->k_/this->epsilon_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
buoyantkEpsilonRealizableBCFDK<BasicTurbulenceModel>
    ::buoyantkEpsilonRealizableBCFDK
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& type
)
:
    buoyantkEpsilonBaseBCFDK
    <
        realizableKE,
        tmp<volScalarField>,
        BasicTurbulenceModel
    >
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        type
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
