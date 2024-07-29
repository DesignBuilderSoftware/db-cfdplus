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

#include "buoyantkEpsilonBaseBCFDK.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template
<
    template<class> class kEpsilonBaseModel,
    class cmuType,
    class BasicTurbulenceModel
>
buoyantkEpsilonBaseBCFDK<kEpsilonBaseModel,cmuType,BasicTurbulenceModel>
    ::buoyantkEpsilonBaseBCFDK
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
    kEpsilonBaseModel<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        type
    ),
    Cbe_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cbe",
            this->coeffDict_,
            1.0
        )
    ),
    Cbk_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cbk",
            this->coeffDict_,
            1.0
        )
    ),
    Tanhe_
    (
        Switch::lookupOrAddToDict
        (
            "Tanhe",
            this->coeffDict_,
            false
        )
    ),
    Tanhk_
    (
        Switch::lookupOrAddToDict
        (
            "Tanhk",
            this->coeffDict_,
            false
        )
    )
{
    /* Not implemented here
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
    */
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template
<
    template<class> class kEpsilonBaseModel,
    class cmuType,
    class BasicTurbulenceModel
>
bool buoyantkEpsilonBaseBCFDK<kEpsilonBaseModel,cmuType,BasicTurbulenceModel>
    ::read()
{
    if (kEpsilonBaseModel<BasicTurbulenceModel>::read())
    {
        Cbe_.readIfPresent(this->coeffDict());
        Cbk_.readIfPresent(this->coeffDict());
        Tanhe_.readIfPresent("Tanhe", this->coeffDict());
        Tanhk_.readIfPresent("Tanhk", this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template
<
    template<class> class kEpsilonBaseModel,
    class cmuType,
    class BasicTurbulenceModel
>
tmp<volScalarField>
buoyantkEpsilonBaseBCFDK<kEpsilonBaseModel,cmuType,BasicTurbulenceModel>
    :: tanhU
(
    const uniformDimensionedVectorField& g
) const
{
    vector gHat(g.value()/mag(g.value()));

    volScalarField v(gHat & this->U_);
    volScalarField u
    (
        mag(this->U_ - gHat*v)
      + dimensionedScalar(dimVelocity, small)
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "tanhke",
            tanh(mag(v)/u)
        )
    );
}


template
<
    template<class> class kEpsilonBaseModel,
    class cmuType,
    class BasicTurbulenceModel
>
tmp<volScalarField>
buoyantkEpsilonBaseBCFDK<kEpsilonBaseModel,cmuType,BasicTurbulenceModel>
    ::Gcoef(const uniformDimensionedVectorField& g) const
{
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

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "Pb",
            mag
            (
                this->rho_ * Cmu() * this->k_ * beta
              /
                Prt * (g & fvc::grad(T))
            )
          /
            (this->epsilon_ + this->epsilonMin_)
        )
    );
}

template
<
    template<class> class kEpsilonBaseModel,
    class cmuType,
    class BasicTurbulenceModel
>
tmp<fvScalarMatrix>
buoyantkEpsilonBaseBCFDK<kEpsilonBaseModel,cmuType,BasicTurbulenceModel>
    ::kSource() const
{
    const uniformDimensionedVectorField& g =
        this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");

    if (mag(g.value()) > small)
    {
        tmp<volScalarField> Gb(Cbk_*Gcoef(g));

        if(Tanhk_)
        {
            Gb.ref() *= tanhU(g);
        }

        return fvm::SuSp(Gb, this->k_);
    }
    else
    {
        return kEpsilonBaseModel<BasicTurbulenceModel>::kSource();
    }
}


template
<
    template<class> class kEpsilonBaseModel,
    class cmuType,
    class BasicTurbulenceModel
>
tmp<fvScalarMatrix>
buoyantkEpsilonBaseBCFDK<kEpsilonBaseModel,cmuType,BasicTurbulenceModel>
    ::epsilonSource() const
{
    const uniformDimensionedVectorField& g =
        this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");

    if (mag(g.value()) > small)
    {
        tmp<volScalarField> Gb(Cbe_*Gcoef(g));

        if(Tanhe_)
        {
            Gb.ref() *= tanhU(g);
        }

        //WARNING: C1_ is not used here due to the original request from the
        //  client, which is by we aren't multiplying by C1_
        return fvm::SuSp(Gb,  this->epsilon_);
    }
    else
    {
        return kEpsilonBaseModel<BasicTurbulenceModel>::epsilonSource();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
