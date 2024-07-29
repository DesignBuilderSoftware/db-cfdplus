/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's ManualInjection class.
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

#include "mathematicalConstants.H"
#include "ParticlesAsRaysInjection.T.H"
#include "PackedBoolList.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticlesAsRaysInjection<CloudType>::ParticlesAsRaysInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
    particleDiameter_
    (
        readScalar(this->coeffDict().lookup("particleDiameter"))
    ),
    particlePositionsFieldName_
    (
        this->coeffDict().lookup("positionFieldName")
    ),
    particleVelocitiesFieldName_
    (
        this->coeffDict().lookup("velocityFieldName")
    ),
    particlePositions_(),
    particleVelocities_(),
    injectorCells_(),
    injectorTetFaces_(),
    injectorTetPts_()
{
    updateMesh();

    // Determine volume of particles to inject
    this->volumeTotal_ =
        pow3(particleDiameter_)*pi/6.0
      * injectorCells_.size();
}


template<class CloudType>
Foam::ParticlesAsRaysInjection<CloudType>::ParticlesAsRaysInjection
(
    const ParticlesAsRaysInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    particleDiameter_(im.particleDiameter_),
    particlePositionsFieldName_(im.particlePositionsFieldName_),
    particleVelocitiesFieldName_(im.particleVelocitiesFieldName_),
    particlePositions_(im.particlePositions_),
    particleVelocities_(im.particleVelocities_),
    injectorCells_(im.injectorCells_),
    injectorTetFaces_(im.injectorTetFaces_),
    injectorTetPts_(im.injectorTetPts_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticlesAsRaysInjection<CloudType>::~ParticlesAsRaysInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticlesAsRaysInjection<CloudType>::updateMesh()
{
    //- Initial particle positions field
    const volVectorField &initialParticlePositions_
    (
        this->owner().db().objectRegistry::template
        lookupObject<volVectorField>
        (
            particlePositionsFieldName_
        )
    );

    //- Initial particle velocities field
    const volVectorField &initialParticleVelocities_
    (
        this->owner().db().objectRegistry::template
        lookupObject<volVectorField>
        (
            particleVelocitiesFieldName_
        )
    );

    particlePositions_ = initialParticlePositions_.primitiveField();
    particleVelocities_ = initialParticleVelocities_.primitiveField();

    injectorCells_.resize(particlePositions_.size(), -1);
    injectorTetFaces_.resize(particlePositions_.size(), -1);
    injectorTetPts_.resize(particlePositions_.size(), -1);

    label nRejected = 0;
    PackedBoolList keep(particlePositions_.size(), true);

    // NOTE: Do not use the search algorithm as used bu ManualInjection,
    // because that code is designed to use the same list on all processors,
    // while ours is already split by processor
    forAll(particlePositions_, pI)
    {
        vector pos = particlePositions_[pI];

        this->owner().mesh().findCellFacePt
        (
            pos,
            injectorCells_[pI],
            injectorTetFaces_[pI],
            injectorTetPts_[pI]
        );

        if (injectorCells_[pI] < 0)
        {
            keep[pI] = false;
            nRejected++;
        }
    }

    if (nRejected > 0)
    {
        inplaceSubset(keep, particlePositions_);
        inplaceSubset(keep, particleVelocities_);
        inplaceSubset(keep, injectorCells_);
        inplaceSubset(keep, injectorTetFaces_);
        inplaceSubset(keep, injectorTetPts_);
    }

    //Gather the total number of rejected points
    nRejected = returnReduce(nRejected, sumOp<label>());
    if (nRejected > 0)
    {
        Info<< "WARNING (ParticlesAsRaysInjection):"
            << " Several particles ignored due to flaws in the mesh."
            << " Total number: " << nRejected
            << endl;
    }
}


template<class CloudType>
Foam::scalar Foam::ParticlesAsRaysInjection<CloudType>::timeEnd() const
{
    // Injection is instantaneous - but allow for a finite interval to
    // avoid numerical issues when interval is zero
    return ROOTVSMALL;
}


template<class CloudType>
Foam::label Foam::ParticlesAsRaysInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    if ((0.0 >= time0) && (0.0 < time1))
    {
        return particlePositions_.size();
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::ParticlesAsRaysInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
)
{
    // All parcels introduced at SOI
    if ((0.0 >= time0) && (0.0 < time1))
    {
        return this->volumeTotal_;
    }
    else
    {
        return 0.0;
    }
}


template<class CloudType>
void Foam::ParticlesAsRaysInjection<CloudType>::setPositionAndCell
(
    const label parcelI,
    const label,
    const scalar,
    vector& position,
    label& cellOwner,
    label& tetFacei,
    label& tetPti
)
{
    position = particlePositions_[parcelI];
    cellOwner = injectorCells_[parcelI];
    tetFacei = injectorTetFaces_[parcelI];
    tetPti = injectorTetPts_[parcelI];
}


template<class CloudType>
void Foam::ParticlesAsRaysInjection<CloudType>::setProperties
(
    const label parcelI,
    const label,
    const scalar,
    typename CloudType::parcelType& parcel
)
{
    // set particle velocity
    parcel.U() = particleVelocities_[parcelI];

    // set particle diameter
    parcel.d() = particleDiameter_;
}


template<class CloudType>
bool Foam::ParticlesAsRaysInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::ParticlesAsRaysInjection<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //
