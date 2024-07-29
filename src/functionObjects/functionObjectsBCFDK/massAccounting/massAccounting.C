/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Derived from OpenFOAM's scalarTransport function object.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2013-2018 FSD blueCAPE Lda  http://www.bluecape.com.pt
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

#include "massAccounting.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

//#include "fanRedirectedOutletVelocityFvPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(massAccounting, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        massAccounting,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::massAccounting::report()
{
    scalar inletPhi = 0.0;
    scalar fixedMassOut = 0.0;

    surfaceScalarField& phi =
        mesh_.lookupObjectRef<surfaceScalarField>(phiName_);

    surfaceScalarField::Boundary& bphi =
        phi.boundaryFieldRef();

    /*
    // Flow rate correction, which for now we won't leave it in, since it
    // didn't work as intended
    volVectorField& U = mesh_.lookupObjectRef<volVectorField>("U");

    volVectorField::Boundary& bU =
        U.boundaryFieldRef();

    forAll(bU, patchi)
    {
        fvPatchVectorField& Up = bU[patchi];
        fvsPatchScalarField& phip = bphi[patchi];

        if (!phip.coupled())
        {
            if (isA<fanRedirectedOutletVelocityFvPatchField>(Up))
            {
                Up.evaluate(Pstream::commsTypes::blocking);

                phip = Up & Up.patch().Sf();

                if (mesh_.foundObject<volScalarField>("rho"))
                {
                    const fvPatchField<scalar>& rhop =
                        Up.patch().lookupPatchField<volScalarField, scalar>("rho");

                    phip *= rhop;
                }

            }
        }
    }
    */

    forAll(bphi, patchi)
    {
        const fvsPatchScalarField& phip = bphi[patchi];

        if (!phip.coupled())
        {
            forAll(phip, pFacei)
            {
                if(phip[pFacei] < 0.0)
                {
                  inletPhi -= phip[pFacei];
                }
                else
                {
                  fixedMassOut += phip[pFacei];
                }
            }
        }
    }

    reduce(inletPhi, sumOp<scalar>());
    reduce(fixedMassOut, sumOp<scalar>());

    scalar imbalance = mag(fixedMassOut - inletPhi);

    // Just to get a notion in debugging on what's going on within the domain
    scalar totalFlux = vSmall + sum(mag(phi)).value();

    Info<< "Inlet phi: " << inletPhi << endl;
    Info<< "Mass flow imbalance: " << imbalance << endl;
    Info<< "Total flux: " << totalFlux << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::massAccounting::massAccounting
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::massAccounting::~massAccounting()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::massAccounting::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    phiName_ = dict.lookupOrDefault<word>("phi", "phi");

    return true;
}


bool Foam::functionObjects::massAccounting::execute()
{
    return true;
}


bool Foam::functionObjects::massAccounting::write()
{
    report();

    return true;
}


// ************************************************************************* //
