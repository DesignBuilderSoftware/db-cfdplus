Info<< "Constructing kinematicCloud " << kinematicCloudName+name(i) << endl;

// Notes:
// - When running in parallel, there is an independent particle ID counter per 
//   processor, where each particle's identifier is retrievable with origId()
//   and fortunately there is an auxiliary ID for processor origProc()
// - The counter is this static public variable: Foam::particle::particleCount_
// - We reset this counter when going through here, otherwise it's not possible
//   to have a one-to-one relation between the original particle location and
//   the cell centre.

Foam::particle::particleCount_ = 0;

basicKinematicCollidingCloud kinematicCloud
(
    kinematicCloudName,
    rhoInf,
    U,
    mu,
    g
);

kinematicCloud.rename(kinematicCloudName+name(i));

