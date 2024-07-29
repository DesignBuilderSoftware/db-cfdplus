
{
    label counter = 0;
    forAllIter(basicKinematicCollidingCloud, kinematicCloud, pIter)
    {
        basicKinematicCollidingCloud::parcelType& p = pIter();

        if (p.active())
        {
            counter++;
        }
    }

    counter = returnReduce(counter, sumOp<label>());

    if (counter == 0)
    {
        Info<< "All particles have collided." << endl;
        break;
    }
}
