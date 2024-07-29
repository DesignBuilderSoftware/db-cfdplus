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

Class
    Foam::ParticlesAsRaysInjection

Description
    Predefined injection, given through existing position and velocity fields.
    Diameter is fixed.

    Example:
    \verbatim
        particlesAsRays
        {
            type            particlesAsRays;
            massTotal       0;
            parcelBasisType fixed;
            nParticle       1;
            SOI             0;
            particleDiameter      0.00001;
            positionFieldName   "particlePositions";
            velocityFieldName   "particleVelocities";
        }
    \endverbatim


SourceFiles
    ParticlesAsRaysInjection.C

\*---------------------------------------------------------------------------*/

#ifndef ParticlesAsRaysInjection_H
#define ParticlesAsRaysInjection_H

#if defined(WIN32) || defined(WIN64)
#include "InjectionModel.T.H"
#else
#include "InjectionModel.H"
#endif

#include "Switch.H"
#include "vectorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ParticlesAsRaysInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParticlesAsRaysInjection
:
    public InjectionModel<CloudType>
{
    // Private Data

        //- Fixed diameter for all particle [m]
        const scalar particleDiameter_;

        //- Name of field that has the initial particle positions
        const word particlePositionsFieldName_;

        //- Name of field that has the initial particle velocities
        const word particleVelocitiesFieldName_;

        //- Initial particle positions field
        vectorList particlePositions_;

        //- Initial particle velocities field
        vectorList particleVelocities_;

        //- List of cell labels corresponding to injector positions
        labelList injectorCells_;

        //- List of tetFace labels corresponding to injector positions
        labelList injectorTetFaces_;

        //- List of tetPt labels corresponding to injector positions
        labelList injectorTetPts_;


public:

    //- Runtime type information
    TypeName("particlesAsRays");


    // Constructors

        //- Construct from dictionary
        ParticlesAsRaysInjection
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        ParticlesAsRaysInjection(const ParticlesAsRaysInjection<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType>> clone() const
        {
            return autoPtr<InjectionModel<CloudType>>
            (
                new ParticlesAsRaysInjection<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ParticlesAsRaysInjection();


    // Member Functions

        //- Set injector locations when mesh is updated
        virtual void updateMesh();

        //- Return the end-of-injection time
        scalar timeEnd() const;

        //- Number of parcels to introduce relative to SOI
        virtual label parcelsToInject(const scalar time0, const scalar time1);

        //- Volume of parcels to introduce relative to SOI
        virtual scalar volumeToInject(const scalar time0, const scalar time1);


        // Injection geometry

            //- Set the injection position and owner cell, tetFace and tetPt
            virtual void setPositionAndCell
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                vector& position,
                label& cellOwner,
                label& tetFacei,
                label& tetPti
            );

            //- Set the parcel properties
            virtual void setProperties
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                typename CloudType::parcelType& parcel
            );

            //- Flag to identify whether model fully describes the parcel
            virtual bool fullyDescribed() const;

            //- Return flag to identify whether or not injection of parcelI is
            //  permitted
            virtual bool validInjection(const label parcelI);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ParticlesAsRaysInjection.T.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
