/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's moveDynamicMesh.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2014-2016 FSD blueCAPE Lda  http://www.bluecape.com.pt
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
    applyDynamicMesh

Description
    Mesh motion and topological mesh changes utility that writes to the place
    where the mesh came from.

\*---------------------------------------------------------------------------*/
#if defined( WIN32 ) || defined( WIN64 )
#include "Time.T.H"
#else
#include "Time.H"
#endif

#include "argList.H"
#include "fvMesh.H"
#include "dynamicFvMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addOverwriteOption.H"
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedDynamicFvMesh.H"

    const bool overwrite = args.optionFound("overwrite");
    const word oldInstance = mesh.pointsInstance();

    Info<< "Time = " << runTime.timeName() << endl;

    //Force the mesh motion to occur
    mesh.update();

    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }

    //Take out the dynamic motion items, without affecting the mesh itself
    mesh.clearOut();

    //save the mesh
    mesh.write();

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
