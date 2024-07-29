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

#include "scalarTransportBCFDK.H"
#include "surfaceFields.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "kinematicMomentumTransportModel.H"
#include "fluidThermoMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(scalarTransportBCFDK, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        scalarTransportBCFDK,
        dictionary
    );

    template<>
    const char* NamedEnum
    <
        scalarTransportBCFDK::diffusionTypes,
        3
    >::names[] =
    {
        "constant",
        "existing",
        "inferred"
    };
}
}

const Foam::NamedEnum
<
    Foam::functionObjects::scalarTransportBCFDK::diffusionTypes,
    3
> Foam::functionObjects::scalarTransportBCFDK::diffusionTypeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::functionObjects::scalarTransportBCFDK::D
(
    const surfaceScalarField& phi
) const
{
    typedef incompressible::momentumTransportModel icoModel;
    typedef compressible::momentumTransportModel cmpModel;

    word Dname("DT_" + fieldName_);

    if (dType_ == constant)
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    Dname,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("localD", phi.dimensions()/dimLength, D_)
            )
        );
    }
    else if(dType_ == existing)
    {
        return mesh_.lookupObject<volScalarField>(Dname);
    }
    else if (mesh_.foundObject<icoModel>(momentumTransportModel::typeName))
    {
        const icoModel& model = mesh_.lookupObject<icoModel>
        (
            momentumTransportModel::typeName
        );

        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    Dname,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                alphaD_*model.nu() + alphaDt_*model.nut()
            )
        );
    }
    else if (mesh_.foundObject<cmpModel>(momentumTransportModel::typeName))
    {
        const cmpModel& model = mesh_.lookupObject<cmpModel>
        (
            momentumTransportModel::typeName
        );

        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    Dname,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                alphaD_*model.mu() + alphaDt_*model.mut()
            )
        );
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    Dname,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar
                (
                    "alphaDzero",
                    phi.dimensions()/dimLength,
                    0.0
                )
            )
        );
    }
}


bool Foam::functionObjects::scalarTransportBCFDK::stoppingCriteria
(
    label iter,
    scalar& resd,
    scalar& resd_0,
    const scalar initialResidual
)
{
    scalar rdrop = 1.0;
    scalar normFactor = 1.0;

    {
        const dictionary& matrixSolverDict = mesh_.solverDict(schemesField_);
        word normFactorName = "normFactor";

        if (iter == 0)
        {
            if (matrixSolverDict.found(normFactorName))
            {
                matrixSolverDict.lookup(normFactorName) >> normFactor;

                if (debug)
                {
                    Info<< "Normalization factor for " << fieldName_
                        << " is " << normFactor
                        << endl;
                }
            }
        }
        else if (iter > 0)
        {
            const word normFactorNameList(normFactorName + "List");
            scalarList extraNormFactors;

            if (matrixSolverDict.found(normFactorNameList))
            {
                matrixSolverDict.lookup(normFactorNameList) >>
                    extraNormFactors;

                if((iter-2) < extraNormFactors.size())
                {
                    normFactor = extraNormFactors[iter-1];
                }
            }
        }
    }

    // Initialize the residual monitoring
    if(iter == 0)
    {
        resd = initialResidual*normFactor;
        resd_0 = resd;
    }
    else
    {
        resd = initialResidual*normFactor;
        rdrop = resd/(resd_0 + VSMALL);
    }

    return rdrop <= residualDropTol_ || resd <= minimumResidual_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::scalarTransportBCFDK::scalarTransportBCFDK
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldName_(dict.lookupOrDefault<word>("field", "s")),
    D_(0),
    dType_(diffusionTypeNames_.read(dict.lookup("diffusionType"))),
    nCorr_(0),
    fvOptions_(mesh_),
    createNewScalar_(dict.lookupOrDefault<Switch>("newField", true))
{
    if (createNewScalar_)
    {
        store(
            tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        fieldName_,
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_
                )
            )
        );
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::scalarTransportBCFDK::~scalarTransportBCFDK()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::scalarTransportBCFDK::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    debug = dict.lookupOrDefault<label>("debug", debug);

    phiName_ = dict.lookupOrDefault<word>("phi", "phi");
    rhoName_ = dict.lookupOrDefault<word>("rho", "rho");
    schemesField_ = dict.lookupOrDefault<word>("schemesField", fieldName_);

    alphaD_ = dict.lookupOrDefault("alphaD", 1.0);
    alphaDt_ = dict.lookupOrDefault("alphaDt", 1.0);

    minFieldValueIsSet_ = dict.readIfPresent("minValue", minFieldValue_);
    resetFieldValueIsSet_ = dict.readIfPresent("resetValue", resetFieldValue_);

    dict.readIfPresent("nCorr", nCorr_);
    residualDropTol_ = dict.lookupOrDefault<scalar>("residualDropTol", 1.E-4);
    minimumResidual_ = dict.lookupOrDefault<scalar>("minimumResidual", 1.E-6);

    if (dict.found("fvOptions"))
    {
        fvOptions_.reset(dict.subDict("fvOptions"));
    }

    return true;
}


bool Foam::functionObjects::scalarTransportBCFDK::execute()
{
    Info<< type() << " " << name() << " execute:" << endl;

    const surfaceScalarField& phi =
        mesh_.lookupObject<surfaceScalarField>(phiName_);

    // Calculate the diffusivity
    Foam::tmp<Foam::volScalarField> tD(this->D(phi));
    const volScalarField &D = tD();

    // Access the scalar field to operate on
    volScalarField& s_ = mesh_.lookupObjectRef<volScalarField>(fieldName_);

    // Reset the field, if so requested
    if(resetFieldValueIsSet_)
    {
        s_ == dimensionedScalar
        (
            "zero_"+fieldName_,
            s_.dimensions(),
            resetFieldValue_
        );
    }

    word divScheme("div(phi," + schemesField_ + ")");
    word laplacianScheme("laplacian(" + D.name() + "," + schemesField_ + ")");

    // Set under-relaxation coeffs
    scalar relaxCoeff = 0.0;
    if (mesh_.relaxEquation(schemesField_))
    {
        relaxCoeff = mesh_.equationRelaxationFactor(schemesField_);
    }

    scalar relaxFieldCoeff = 0.0;
    if (mesh_.relaxField(schemesField_))
    {
        relaxFieldCoeff = mesh_.fieldRelaxationFactor(schemesField_);
    }

    scalar resd = -1;
    scalar resd_0 = -1;
    solverPerformance performance;

    if (phi.dimensions() == dimMass/dimTime)
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);

        for (int i = 0; i <= nCorr_; i++)
        {
            s_.storePrevIter();

            fvScalarMatrix sEqn
            (
                fvm::ddt(rho, s_)
              + fvm::div(phi, s_, divScheme)
              - fvm::laplacian(D, s_, laplacianScheme)
             ==
                fvOptions_(rho, s_)
            );

            sEqn.relax(relaxCoeff);

            fvOptions_.constrain(sEqn);

            performance = sEqn.solve(mesh_.solverDict(schemesField_));

            fvOptions_.correct(s_);

            // Explicit relaxation
            s_.relax(relaxFieldCoeff);

            if(debug)
            {
                volScalarField tmpS
                (
                    IOobject
                    (
                        fieldName_ + "_iteration" + Foam::name(i),
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    s_
                );
                tmpS.write();
            }

            if
            (
                stoppingCriteria
                (
                    i,
                    resd,
                    resd_0,
                    performance.initialResidual()
                )
            ) break;
        }
    }
    else if (phi.dimensions() == dimVolume/dimTime)
    {
        for (int i = 0; i <= nCorr_; i++)
        {
            s_.storePrevIter();

            fvScalarMatrix sEqn
            (
                fvm::ddt(s_)
              + fvm::div(phi, s_, divScheme)
              - fvm::laplacian(D, s_, laplacianScheme)
             ==
                fvOptions_(s_)
            );

            sEqn.relax(relaxCoeff);

            fvOptions_.constrain(sEqn);

            performance = sEqn.solve(mesh_.solverDict(schemesField_));

            fvOptions_.correct(s_);

            // Explicit relaxation
            s_.relax(relaxFieldCoeff);

            if(debug)
            {
                volScalarField tmpS
                (
                    IOobject
                    (
                        fieldName_ + "_iteration" + Foam::name(i),
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    s_
                );
                tmpS.write();
            }

            if
            (
                stoppingCriteria
                (
                    i,
                    resd,
                    resd_0,
                    performance.initialResidual()
                )
            ) break;
        }
    }
    else
    {
        FatalErrorInFunction
            << "Incompatible dimensions for phi: " << phi.dimensions() << nl
            << "Dimensions should be " << dimMass/dimTime << " or "
            << dimVolume/dimTime << exit(FatalError);
    }

    if(minFieldValueIsSet_)
    {
        s_ = max
        (
            s_,
            dimensionedScalar
            (
                "min_" + fieldName_,
                s_.dimensions(),
                minFieldValue_
            )
        );
    }

    Info<< endl;

    return true;
}


bool Foam::functionObjects::scalarTransportBCFDK::write()
{
    volScalarField& s_ = mesh_.lookupObjectRef<volScalarField>(fieldName_);
    s_.write();

    return true;
}


// ************************************************************************* //
