/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "buoyantkOmegaSSTBCFDK.H"
#include "uniformDimensionedFields.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
buoyantkOmegaSSTBCFDK<BasicTurbulenceModel>::buoyantkOmegaSSTBCFDK
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
    kOmegaSST<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        type
    ),
    Cb_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cb",
            this->coeffDict_,
            1.0
        )
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool buoyantkOmegaSSTBCFDK<BasicTurbulenceModel>::read()
{
    if (kOmegaSST<BasicTurbulenceModel>::read())
    {
        Cb_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
tmp<volScalarField>
buoyantkOmegaSSTBCFDK<BasicTurbulenceModel>::Gcoef() const
{
    const uniformDimensionedVectorField& g =
        this->mesh_.objectRegistry::
        lookupObject<uniformDimensionedVectorField>("g");

    const IOdictionary& transportProperties =
        this->mesh_.objectRegistry::
        lookupObject<IOdictionary>("transportProperties");

    const volScalarField& T =
        this->mesh_.objectRegistry::
        lookupObject<volScalarField>("T");

    const dimensionedScalar Prt("Prt", dimless, transportProperties);
    const dimensionedScalar beta
    (
        "beta",
        dimless/dimTemperature,
        transportProperties
    );

    // Note: The old implementation used "nut_", but we went with the
    // convention from OpenFOAM 5 instead, namely that:
    //    this->nut_ = Cmu_*sqr(k_)/epsilon_;
    // where as in k-omega SST, it's this:
    //    nut_ = a1_*k_/max(a1_*omega_, b1_*F23()*sqrt(S2));

    tmp<volTensorField> tgradU = fvc::grad(this->U_);
    volScalarField S2(2*magSqr(symm(tgradU())));
    tgradU.clear();

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "Pb",
            mag(this->rho_ * this->a1_ * beta / Prt * (g & fvc::grad(T)))
          /
            max(this->a1_*this->omega_, this->b1_*this->F23()*sqrt(S2))
        )
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix>
buoyantkOmegaSSTBCFDK<BasicTurbulenceModel>::kSource() const
{
    const uniformDimensionedVectorField& g =
        this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");

    if (mag(g.value()) > small)
    {
        return fvm::SuSp(Gcoef(), this->k_);
    }
    else
    {
        return kOmegaSST<BasicTurbulenceModel>::kSource();
    }
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix>
buoyantkOmegaSSTBCFDK<BasicTurbulenceModel>::omegaSource() const
{
    const uniformDimensionedVectorField& g =
        this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");

    if (mag(g.value()) > small)
    {
        const volScalarField CDkOmega
        (
            (2*this->alphaOmega2_)
          * (
                fvc::grad(this->k_) & fvc::grad(this->omega_)
            )/this->omega_
        );

        const volScalarField F1(this->F1(CDkOmega));

        //Original implementation we had in 2.2.x was:
        //omega_/k_*(scalar(5.0/9.0)*F1+(scalar(1)-F1)*scalar(4.0/9.0))*Cb_*Pb
        //
        //Which according to the notes in the task BCFDK-303, where the large
        //parenthesis is part of the SST mixing mechanism.

        return fvm::SuSp
        (
            (scalar(5.0/9.0)*F1+(scalar(1)-F1)*scalar(4.0/9.0))*Cb_*Gcoef(),
            this->omega_
        );
    }
    else
    {
        return kOmegaSST<BasicTurbulenceModel>::omegaSource();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
