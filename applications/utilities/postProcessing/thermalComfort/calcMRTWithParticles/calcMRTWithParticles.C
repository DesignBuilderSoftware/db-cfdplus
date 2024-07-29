/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's icoUncoupledKinematicParcelFoam.
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

Application
    calcMRTWithParticles

Description
    Utility for calculating the Mean Radiant Temperature field, through the use
    of transient passive transport of kinematic particle clouds for ray
    tracing, by ignoring the flow field and fluid properties.

    Uses a pre-calculated velocity field to evolve the cloud.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "basicKinematicCollidingCloud.H"
#include "makeParticlesAsRaysInjectionModel.H"
#include "meshSearch.H"
#include "controlledTime/controlledTime.H"
#include "parallelPatchCloud/parallelPatchCloudSet.H"

#if defined(WIN32) || defined(WIN64)
#include "Tuple2.T.H"
#include "List.T.H"
#else
#include "Tuple2.H"
#include "List.H"
#endif

#include "mapDistribute.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*
    Note: the points must be entered in an order that the lines that connects
    them do not cross. For instance:
    D ______ C
    |       |
    |       |
    |_______|
    B       A
    the order, in the calling function, can not be A, B, C and D but yes A, B,
    D and C.
*/


scalar calcAreaFrom4Points
(
    const point & A,
    const point & B,
    const point & C,
    const point & D,
    const point & centroid
)
{
    scalar area;
    scalar sumA = 0.0;

    const label nPoints = 4;
    List<point> p(nPoints);
    p[0] = A;
    p[1] = B;
    p[2] = C;
    p[3] = D;

    for (label i = 0; i < nPoints; i++)
    {
        const point& nextPoint = p[(i + 1) % nPoints];

        vector n = (nextPoint - p[i])^(centroid - p[i]);
        scalar a = mag(n);

        sumA += a;
    }

    // This is to deal with zero-area faces. Mark very small faces
    // to be detected in e.g., processorPolyPatch.
    if (sumA < ROOTVSMALL)
    {
        area = 0.0;
    }
    else
    {
        area = 0.5*sumA;
    }

    return area;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "cloudName",
        "name",
        "specify alternative cloud name. default is 'kinematicCloud'"
    );

    argList::addBoolOption
    (
        "debug",
        "Debug mode"
    );

    timeSelector::addOptions();

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const bool debug = args.optionFound("debug");
    const word dictName("calcMRTWithParticlesDict");

    // Reading calcMRTWithParticlesDict
    IOdictionary calcMRTWithParticlesDict
    (
       IOobject
       (
            "calcMRTWithParticlesDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
       )
    );

    #include "createFields.H"
    #include "createTimingVariables.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    const label numberOfSegments = 
        calcMRTWithParticlesDict.lookup<label>("numberOfSegments");
    const label numberOfStripes = numberOfSegments/2;
    const label numberOfPatches = numberOfSegments*numberOfStripes;

    // Calculating the area of each patch and direction of each particle
    #include "calcAreasAndCentroids.H"

    runTime.functionObjects().start();

    for(label i=0; i < numberOfPatches; i++)
    {
        Info
            << nl << nl
            << "Patch number " << i << " out of " << numberOfPatches
            << nl << endl;

        // Define the particle directions (as a velocity vector)
        particleVelocities = characteristicLength*patchCentroids[i];

        // To be used for ensuring the calculations are done right
        #include "resetAndUpdateTimings.H"

        // To be used for each ray direction
        #include "createParticleField.H"


        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        Info<< "\nStarting time loop\n" << endl;

        while (runTime.loop())
        {
            Info<< "Time = " << runTime.timeName() << nl << endl;

            Info<< "Evolving " << kinematicCloud.name() << endl;

            kinematicCloud.evolve();

            if (debug)
            {
                runTime.write();
            }

            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;

            #include "particleCloudStoppingCriteria.H"

        }

        if (debug)
        {
            runTime.writeNow();
        }

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // Sample this cloud on all patches
        #include "sampleParticlesOnPatches.H"

        // Transfer sampled particle result back to their original processors
        #include "transferSamplesBackToOriginalProcessors.H"

        // MRT calculation
        #include "addIncrementalMRT.H"

    }

    runTime.setTime(originalTimeSnapshot, originalTimeSnapshotIndex);
    sphereAreas.boundaryFieldRef() = 1.0;
    MRT = MRT/(sphereAreas+SMALL);
    MRT.write();

    runTime.functionObjects().execute();
    runTime.functionObjects().end();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
