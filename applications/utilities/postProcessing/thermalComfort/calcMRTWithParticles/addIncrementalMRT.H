{
    volScalarField localPatchAreas(
        patchAreas[i]*pow(characteristicLength, 2.0)
        );

    MRT += MRTp*localPatchAreas;

    // Do not unaccount unused areas, because there are situations where some
    // cells are not able to have a single sample
    sphereAreas += neg(-MRTp)*localPatchAreas;

    if (debug)
    {
        runTime.setTime(originalTimeSnapshot, originalTimeSnapshotIndex);
        MRTp.rename("MRTp"+name(i));
        MRTp.write();        
    }
}

