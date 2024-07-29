Info<< "\nGenerating dummy fields...\n" << endl;

Info<< "Creating dummy field rhoInf\n" << endl;
volScalarField rhoInf
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rho", dimMass/dimVolume, 1.0)
);

Info<< "Creating dummy field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("u", dimless, vector::zero) //dimVelocity
);

Info<< "Creating dummy field mu" << endl;
volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("mu", dimMass/dimLength/dimTime, 1.0)
);

Info<< "Creating dummy vector field g" << endl;
uniformDimensionedVectorField g
(
    IOobject
    (
        "g",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    dimensionedVector("g", dimVelocity/dimTime, vector::zero)
);

// Particle related fields and settings...

word kinematicCloudName("kinematicCloud");
args.optionReadIfPresent("cloudName", kinematicCloudName);

Info<< "Creating field particleVelocities\n" << endl;
volVectorField particleVelocities
(
    IOobject
    (
        "particleVelocities",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    U
);

Info<< "Creating field particlePositions\n" << endl;
volVectorField particlePositions
(
    IOobject
    (
        "particlePositions",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh.C()
);

Info<< "Creating field Mean Radiant Temperature\n" << endl;
volScalarField MRT
(
    IOobject
    (
        "MRT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("MRT", dimTemperature, 0.0)
);

volScalarField MRTp
(
    IOobject
    (
        "MRTp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    MRT
);

volScalarField sphereAreas
(
    IOobject
    (
        "sphereAreas",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("sphereAreas", dimless, 0.0)
);

Info<< "Reading field T\n" << endl;
const volScalarField T
(
    IOobject
    (
        calcMRTWithParticlesDict.lookup<word>("T"),
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volScalarField characteristicLength
(
    IOobject
    (
        "characteristicLength",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("characteristicLength", dimless, 0.0)
);

characteristicLength.primitiveFieldRef() = cbrt(mesh.V());

if (debug)
{
    characteristicLength.write();
}
