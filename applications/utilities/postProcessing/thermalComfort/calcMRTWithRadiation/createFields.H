Info<< "Reading thermophysical properties\n" << endl;

autoPtr<rhoThermo> pThermo(rhoThermo::New(mesh));

rhoThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

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
