// Next step, implement transfering the particles+values back to their
// original positions.
// Note: This is using code adapted from "applications/test/parallel"

// First, we repackage the know values, either for later using in serial mode
// or for swapping between processors.

// origProcID, origID, sampleT
List<Tuple2<label, Tuple2<label, scalar> > > samplingData
(
    kinematicCloud.size()
);

// Package data to be transferred, using the structure above
counter = 0;
forAllIter(basicKinematicCollidingCloud, kinematicCloud, pIter)
{
    basicKinematicCollidingCloud::parcelType& p = pIter();

    samplingData[counter].first() = p.origProc();
    samplingData[counter].second().first() = p.origId();
    samplingData[counter].second().second() = sampledT[counter];

    counter++;
}

if (Pstream::parRun())
{
    // Need to reschuffle the values list back to their original owner
    // positions

    // Count how many to send
    labelList nSend(Pstream::nProcs(), 0);
    forAll(samplingData, i)
    {
        label procI = samplingData[i].first();
        nSend[procI]++;
    }

    // Sync how many to send
    labelListList allNTrans(Pstream::nProcs());
    allNTrans[Pstream::myProcNo()] = nSend;
    combineReduce(allNTrans, UPstream::listEq());

    // Collect items to be sent
    labelListList sendMap(Pstream::nProcs());
    forAll(sendMap, procI)
    {
        sendMap[procI].setSize(nSend[procI]);
    }
    nSend = 0;
    forAll(samplingData, i)
    {
        label procI = samplingData[i].first();
        sendMap[procI][nSend[procI]++] = i;
    }

    // Collect items to be received
    labelListList recvMap(Pstream::nProcs());
    forAll(recvMap, procI)
    {
        recvMap[procI].setSize(allNTrans[procI][Pstream::myProcNo()]);
    }

    label constructSize = 0;
    // Construct with my own elements first
    forAll(recvMap[Pstream::myProcNo()], i)
    {
        recvMap[Pstream::myProcNo()][i] = constructSize++;
    }
    // Construct from other processors
    forAll(recvMap, procI)
    {
        if (procI != Pstream::myProcNo())
        {
            forAll(recvMap[procI], i)
            {
                recvMap[procI][i] = constructSize++;
            }
        }
    }


    // Construct distribute map (destructively)
    mapDistribute map(constructSize, move(sendMap), move(recvMap));

    // Distribute samplingData
    map.distribute(samplingData);
}

// Reset incremental field
MRTp = dimensionedScalar("MRT", dimTemperature, 0.0);

// Place the values on the MRTp field
forAll(samplingData, sampleI)
{
    const label cellI = samplingData[sampleI].second().first();
    const scalar temperatureI = samplingData[sampleI].second().second();

    if(cellI>=0)
    {
      MRTp[cellI] = temperatureI;
    }
}

