/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    yPlusBasedOnWSS

Description
    Calculates and reports y+ for all patches, for thespecified times when 
    using RAS turbulence models.

    Does not work for compressible,

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "incompressible/RAS/RASModel/RASModel.H"

#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void calcIncompressible
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    volVectorField& wallShearStress,
    volScalarField& yPlus
)
{
    #include "createPhi.H"

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::RASModel> model
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );

    const volSymmTensorField Reff(model->devReff());

    forAll(wallShearStress.boundaryField(), patchI)
    {
        wallShearStress.boundaryField()[patchI] =
        (
           -mesh.Sf().boundaryField()[patchI]
           /mesh.magSf().boundaryField()[patchI]
        ) & Reff.boundaryField()[patchI];

        yPlus.boundaryField()[patchI] = 
        model->y()[patchI]
      * sqrt(mag(wallShearStress.boundaryField()[patchI]))
      / model->nu();
    }
}

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "addRegionOption.H"

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();
/*
        Info<< "Calculating wall distance\n" << endl;
        wallDist y(mesh, true);
        Info<< "Writing wall distance to field " << y.name() << nl << endl;
        y.write();
*/

        volScalarField yPlus
        (
            IOobject
            (
                "yPlusWSS",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("yPlus", dimless, 0.0)
        );

        volVectorField wallShearStress
        (
            IOobject
            (
                "wallShearStress",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedVector
            (
                "wallShearStress",
                sqr(dimLength)/sqr(dimTime),
                vector::zero
            )
        );

        IOobject UHeader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (UHeader.headerOk())
        {
            Info<< "Reading field U\n" << endl;
            volVectorField U(UHeader, mesh);

            calcIncompressible(mesh, runTime, U, wallShearStress, yPlus);
        }
        else
        {
            Info<< "    no U field" << endl;
        }

        Info<< "Writing y+ based on the wall shear stress to field " 
            << yPlus.name() << nl << endl;

        yPlus.write();
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
