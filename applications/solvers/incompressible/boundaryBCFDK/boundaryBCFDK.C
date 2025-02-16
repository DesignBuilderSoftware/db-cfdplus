/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's boundaryFoam.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2014 FSD blueCAPE Lda  http://www.bluecape.com.pt
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
    boundaryBCFDK

Description
    Steady-state solver for incompressible, 1D turbulent flow, typically to
    generate boundary layer conditions at an inlet, for use in a simulation.

    Boundary layer code to calculate the U, k and epsilon distributions.
    Used to create inlet boundary conditions for experimental comparisons
    for which U and k have not been measured.
    Turbulence model is runtime selectable.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "wallFvPatch.H"
#include "makeGraph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "interrogateWallPatches.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        fvVectorMatrix divR(turbulence->divDevReff(U));
        divR.source() = flowMask & divR.source();

        fvVectorMatrix UEqn
        (
            divR == gradP
        );

        UEqn.relax();

        UEqn.solve();


        // Correct driving force for a constant volume flow rate
        dimensionedVector UbarStar = flowMask & U.weightedAverage(mesh.V());

        U += (Ubar - UbarStar);
        gradP += (Ubar - UbarStar)/(1.0/UEqn.A())().weightedAverage(mesh.V());


        turbulence->correct();

        Info<< "Uncorrected Ubar = " << (flowDirection & UbarStar.value())
            << ", pressure gradient = " << (flowDirection & gradP.value())
            << endl;

        #include "evaluateNearWall.H"

        if (runTime.outputTime())
        {
            #include "makeGraphs.H"
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
