// FSD blueCAPE Lda 2018/2023: store the normFactor into the solver's dictionary
{
    // Access our dictionary reference in fvSolution
    dictionary& matrixSolverDictOv =
        const_cast<dictionary&>(solverControls_);

    // Define the names of the variables for the current equation and iteration
    word normFactorName = "normFactor";
    word normFactorTimeName("normFactorTime");

    if (cmpt==1)
    {
        normFactorName += word("Y");
        normFactorTimeName += word("Y");
    }
    else if (cmpt==2)
    {
        normFactorName += word("Z");
        normFactorTimeName += word("Z");
    }
    else if (cmpt>0)
    {
        normFactorName += word("UNKNOWN");
        normFactorTimeName += word("UNKNOWN");
    }

    const word normFactorNameList(normFactorName + "List");

    // Figure out if it's a new time step or the same
    const scalar currentTime(matrix().mesh().thisDb().time().value());
    const scalar prevTime = matrixSolverDictOv.lookupOrDefault<scalar>
    (
        normFactorTimeName,
        scalar(-1.0)
    );
    const bool newTime(prevTime<currentTime);

    //write to dict accordingly
    if (newTime)
    {
        matrixSolverDictOv.add(normFactorTimeName, currentTime, true);
        matrixSolverDictOv.add(normFactorName, normFactor, true);
        matrixSolverDictOv.add(normFactorNameList, scalarList(), true);
    }
    else
    {
        scalarList normFactors = matrixSolverDictOv.lookupOrDefault<scalarList>
        (
            normFactorNameList,
            scalarList()
        );
        normFactors.append(normFactor);
        matrixSolverDictOv.add(normFactorNameList, normFactors, true);
    }

    if (debug)
    {
        Info << "Matrix type solver: " << typeName << ", field: " << fieldName_
            << ", norm factor: " << normFactor << endl;
    }
}
// End of hack ------------------------------------------------------------
