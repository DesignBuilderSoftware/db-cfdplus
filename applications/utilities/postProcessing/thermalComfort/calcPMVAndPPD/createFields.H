    Info<< "Reading field U\n" << endl;
    const volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    Info<< "Reading field T\n" << endl;
    const volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    Info<< "Reading field MRT\n" << endl;
    const volScalarField MRT
    (
        IOobject
        (
            "MRT",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volScalarField PPD
    (
        IOobject
        (
            "PPD",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("PPD", dimless, 0.0)
    );

    volScalarField PMV
    (
        IOobject
        (
            "PMV",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("PMV", dimless, 0.0)
    );

    volScalarField comfortIndex
    (
        IOobject
        (
            "comfortIndex",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("comfortIndex", dimless, 0.0)
    );

    volScalarField OpTempCelsius
    (
        IOobject
        (
            "OpTempCelsius",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("OpTempCelsius", dimTemperature, 0.0)
    );

