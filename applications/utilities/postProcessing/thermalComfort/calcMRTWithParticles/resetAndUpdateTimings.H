// Example:
//   - If minimum velocity magnitude is 0.01;
//   - and diagonal length is 1.5m;
//   - then the deltaTimeStep will be 1.5/0.01 = 150s
//   - and the total time will be startTime + deltaTimeStep

// Hard-coded to 3 times the total simulation time to include a generous safety
// margin, although it should always stop at the end of the first time step

const scalar maximumTimeSteps = 3.0;
const scalar pseudoUnlimitedTimeSteps = maximumTimeSteps*10.0;

minimumVelocity = min(mag(particleVelocities.primitiveField()));
minimumVelocity = returnReduce(minimumVelocity, minOp<scalar>());

deltaTimeStep.value() = ceil(diagonalLength / minimumVelocity);

endTime.value() =
    originalTimeSnapshot.value()
  + deltaTimeStep.value()*maximumTimeSteps;
runTime.setEndTime(endTime);

// Resetting starting time stamp
runTime.setTime(originalTimeSnapshot, originalTimeSnapshotIndex);
runTime.setDeltaT(deltaTimeStep);

if (debug)
{
    //save all time steps
    runTime.setWriteInterval(deltaTimeStep.value());

    Info
        << nl
        << "minimumVelocity: " << minimumVelocity << nl
        << endl;
}
else
{
    runTime.setWriteInterval(deltaTimeStep.value()*pseudoUnlimitedTimeSteps);
}
