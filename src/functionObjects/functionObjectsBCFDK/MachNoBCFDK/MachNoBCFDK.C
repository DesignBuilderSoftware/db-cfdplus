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

#include "MachNoBCFDK.H"
#include "volFields.H"
#include "fluidThermo.H"
#include "addToRunTimeSelectionTable.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(MachNoBCFDK, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        MachNoBCFDK,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Copied/adapted from Foam::functionObjects::MachNo::calc()
bool Foam::functionObjects::MachNoBCFDK::compressibleCalc()
{
    if
    (
        foundObject<volVectorField>(fieldName_)
     && foundObject<fluidThermo>(fluidThermo::dictName)
    )
    {
        const fluidThermo& thermo =
            lookupObject<fluidThermo>(fluidThermo::dictName);

        const volVectorField& U = lookupObject<volVectorField>(fieldName_);

        bool isStored = store
        (
            resultName_,
            mag(U)/sqrt(thermo.gamma()*thermo.p()/thermo.rho())
        );

        if(isStored)
        {
            const volScalarField& MachNo =
                lookupObject<volScalarField>(resultName_);

            Info << "Maximum Mach number: " << gMax(MachNo) << endl;
        }

        return true;
    }
    else
    {
        return false;
    }
}


template<class TempType>
bool Foam::functionObjects::MachNoBCFDK::incompressibleCalc
(
  const TempType &temp
)
{
    const volVectorField& U = lookupObject<volVectorField>(fieldName_);

    //the speed of sound formula "Mach = v/c" is taken from
    //http://en.wikipedia.org/wiki/Speed_of_sound#Speed_in_ideal_gases_and_in_air
    //where:
    // - v is the imposed velocity and c=sqrt(gamma*R*T/Mair) is the speed of
    //   sound
    // - gamma = 1.4
    // - R ~= 8.3145 J.mol-1.K-1
    // - Mair = 0.0289645 kg/mol
    //this formula assumes a constant density

    bool isStored = store
    (
        resultName_,
        mag(U)
      /
        sqrt(gamma_ * Foam::constant::physicoChemical::R * temp / airMass_)
    );

    if(isStored)
    {
        const volScalarField& MachNo =
            lookupObject<volScalarField>(resultName_);

        Info << "Maximum Mach number: " << gMax(MachNo) << endl;
    }

    return isStored;
}


bool Foam::functionObjects::MachNoBCFDK::calc()
{
    if
    (
        foundObject<volVectorField>(fieldName_)
    )
    {
        if(isIsothermal_)
        {
            return incompressibleCalc(TRef_);
        }
        else
        {
            return compressibleCalc();
        }
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::MachNoBCFDK::MachNoBCFDK
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "MachNoBCFDK", "U"),
    TRef_("TRef", dimensionSet(0,0,0,1,0), scalar(0.0)),
    gamma_("gamma", dimensionSet(0,0,0,0,0), scalar(1.4)),
    airMass_
    (
        "airMass",
        dimensionSet(1,0,0,0,-1),
        scalar(0.0289645)
    )
{
    resultName_ = "MachNo";
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::MachNoBCFDK::~MachNoBCFDK()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::MachNoBCFDK::read(const dictionary& dict)
{
    Info<< type() << " read dictionary:" << endl;

    fieldExpression::read(dict);

    dict.lookup("isIsothermal") >> isIsothermal_;

    if(isIsothermal_)
    {
      const dictionary& transportProperties =
          obr_.lookupObject<dictionary>("transportProperties");

        Info<< "    Reading reference temperature TRef\n" << endl;
        transportProperties.lookup("TRef") >> TRef_.value();

        if(transportProperties.found("gamma"))
        {
            Info<< "    Reading adiabatic index\n" << endl;
            transportProperties.lookup("gamma") >> gamma_.value();
        }

        if(transportProperties.found("airMass"))
        {
            Info<< "    Reading mean molar mass for dry air\n" << endl;
            transportProperties.lookup("airMass") >> airMass_.value();
        }
    }

    return true;
}

// ************************************************************************* //
