    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<rhoThermo> pThermo(rhoThermo::New(mesh));

    rhoThermo& thermo = pThermo();
    thermo.validate(args.executable(), "h", "e");
