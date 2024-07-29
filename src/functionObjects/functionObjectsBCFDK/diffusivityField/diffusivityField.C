/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Derived from OpenFOAM's MachNo function object.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2013-2020 FSD blueCAPE Lda  http://www.bluecape.com.pt
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

#include "diffusivityField.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(diffusivityField, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        diffusivityField,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::diffusivityField::calc()
{
    const surfaceScalarField& phi =
        mesh_.lookupObject<surfaceScalarField>(phiName_);

    //This is a bit of a hack for what we want to do, at least for the near
    //future. Generalizing this would require too much work for now
    if(mesh_.foundObject<volScalarField>("nut"))
    {
        const volScalarField& nut =
            mesh_.lookupObject<volScalarField>("nut");

        if (phi.dimensions() == dimMass/dimTime)
        {
            const volScalarField& rho =
                mesh_.lookupObject<volScalarField>(rhoName_);

            return store
            (
                resultName_,
                dimensionedScalar
                (
                    "alphaD",
                    phi.dimensions()/dimLength,
                    alphaD_
                )
                + alphaDt_*nut*rho
            );
        }
        else //if (phi.dimensions() == dimVolume/dimTime)
        {
            return store
            (
                resultName_,
                dimensionedScalar
                (
                    resultName_,
                    phi.dimensions()/dimLength,
                    alphaD_
                )
                + alphaDt_*nut
            );
        }
    }
    else
    {
        WarningInFunction
            << "nut field was not found, using only alphaD as-is."
            << endl;

        return store
        (
            tmp<volScalarField>(new volScalarField
            (
                IOobject
                (
                    resultName_,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar
                (
                    "alphaD",
                    phi.dimensions()/dimLength,
                    alphaD_
                )
            ))
        );
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::diffusivityField::diffusivityField
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "diffusivityField", "DT")
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::diffusivityField::~diffusivityField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::diffusivityField::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    phiName_ = dict.lookupOrDefault<word>("phi", "phi");
    rhoName_ = dict.lookupOrDefault<word>("rho", "rho");
    alphaD_ = dict.lookup<scalar>("alphaD");
    alphaDt_ = dict.lookup<scalar>("alphaDt");

    return true;
}

// ************************************************************************* //
