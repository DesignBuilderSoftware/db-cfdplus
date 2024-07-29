/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's buoyantSimpleFoam.
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
    calcMRTWithRadiation

Description
    Steady-state solver for buoyant, turbulent flow of compressible fluids,
    including radiation, for calculating Mean Radiant Temperature, based on the
    G field.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rhoThermo.H"
#include "radiationModel.H"
#include "simpleControl.H"
#include "wallFvPatch.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        runTime.functionObjects().start();

        Info<< "Time = " << runTime.timeName() << endl;

        mesh.readUpdate();

        #include "createFields.H"
        #include "createRadiationModel.H"

        thermo.correct();
        radiation->correct();

        const volScalarField& G =
            mesh.lookupObject<volScalarField>
            (
                "G"
            );

        MRT = pow(G*0.25/physicoChemical::sigma, 0.25);
        MRT.write();

        runTime.functionObjects().execute();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s" << endl;
    }

    runTime.functionObjects().end();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
