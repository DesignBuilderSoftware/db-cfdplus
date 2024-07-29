Foam::Info<< "Create time\n" << Foam::endl;

Foam::controlledTime runTime(Foam::Time::controlDictName, args);


instantList timeDirs = timeSelector::select0(runTime, args);
if (timeDirs.size() != label(1))
{
    FatalError
        << "createTimingVariables.H: "
        << "There can be only one initial time step." << nl
        << exit(FatalError);
}

const dimensionedScalar originalTimeSnapshot = timeDirs.first().value();
const label originalTimeSnapshotIndex = 
    runTime.findClosestTimeIndex
    (
        runTime.times(),
        originalTimeSnapshot.value()
    );
