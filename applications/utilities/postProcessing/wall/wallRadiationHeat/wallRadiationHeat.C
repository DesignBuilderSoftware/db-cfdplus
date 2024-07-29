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
    wallRadiationHeat

Description
    Calculate heat flux exclusively from radiation on surfaces.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rhoThermo.H"
#include "radiationModel.H"
#include "simpleControl.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    instantList timeDirs = timeSelector::select0(runTime, args);


    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        mesh.readUpdate();

        #include "createFields.H"
        #include "createRadiationModel.H"

        thermo.correct();
        radiation->correct();

        runTime.functionObjects().execute();

        const volScalarField& Qr =
            mesh.lookupObject<volScalarField>
            (
                "qr"
            );

        const volScalarField::Boundary& patchRadHeatFlux =
            Qr.boundaryField();

        const surfaceScalarField::Boundary& magSf =
            mesh.magSf().boundaryField();

        Info<< "\nWall heat fluxes [W]" << endl;
        forAll(patchRadHeatFlux, patchi)
        {
            if (isA<wallFvPatch>(mesh.boundary()[patchi]))
            {
                scalar radFlux = -gSum(magSf[patchi]*patchRadHeatFlux[patchi]);

                Info<< mesh.boundary()[patchi].name() << endl
                    << "    radiative:  " << radFlux << endl;
            }
        }
        Info<< endl;

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s" << endl;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
