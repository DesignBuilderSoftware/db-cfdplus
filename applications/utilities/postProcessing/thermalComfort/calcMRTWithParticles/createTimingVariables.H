// Initialize the time steps
dimensionedScalar endTime = runTime.endTime();
dimensionedScalar deltaTimeStep = runTime.deltaT();

// Calculate chracteristic dimensions
const scalar diagonalLength = boundBox(mesh.points(), true).mag();

if (debug)
{
    Pout
        << nl
        << "diagonalLength: " << diagonalLength << nl
        << endl;
}

// Create minimum velocity variable used for timing calculations
scalar minimumVelocity = 0.0;

// Backup and neutralize the writeInterval
//const scalar originalWriteInterval = runTime.writeInterval();
