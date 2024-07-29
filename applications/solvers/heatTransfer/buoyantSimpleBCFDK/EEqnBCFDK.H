if(thermo.lookupOrAddDefault("disableMechanicalEnergyTerms",false))
{
    Info << "Disabled mechanical energy" << endl;

    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::div(phi, he)
      + thermophysicalTransport->divq(he)
     ==
        radiation->Sh(thermo, he)
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
    radiation->correct();
}
else
{
  #include "EEqn.H"
}