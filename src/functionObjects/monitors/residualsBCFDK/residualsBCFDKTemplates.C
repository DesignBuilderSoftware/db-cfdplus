/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Derived from OpenFOAM's residuals function object.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2023-2023 FSD blueCAPE Lda  http://www.bluecape.com.pt
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "residualsBCFDK.H"
#include "volFields.H"

#if defined( WIN32 ) || defined( WIN64 )
#include "Residuals.T.H"
#include "ListOps.T.H"
#else
#include "Residuals.H"
#include "ListOps.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::residualsBCFDK::writeFileHeader
(
    const label i,
    const word& fieldName,
    const label substep
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (obr_.foundObject<fieldType>(fieldName))
    {
        typename pTraits<Type>::labelType validComponents
        (
            mesh_.validComponents<Type>()
        );

        for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
        {
            if (component(validComponents, cmpt) != -1)
            {
                writeTabbed
                (
                    file(i),
                    fieldName
                  + word(pTraits<Type>::componentNames[cmpt])
                  + subStepName(substep)
                );
            }
        }
    }
}


template<class Type>
void Foam::functionObjects::residualsBCFDK::writeResidual
(
    const word& fieldName,
    const label substep
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (obr_.foundObject<fieldType>(fieldName))
    {
        if (Residuals<Type>::found(mesh_, fieldName))
        {
            const DynamicList<SolverPerformance<Type>>& sp
            (
                Residuals<Type>::field(mesh_, fieldName)
            );

            if (debug)
            {
                Info
                    << "Residuals for " << fieldName << ":" << nl
                    << sp << nl << endl;
            }

            if(substep < sp.size())
            {

                const Type& residual = sp[substep].initialResidual();

                typename pTraits<Type>::labelType validComponents
                (
                    mesh_.validComponents<Type>()
                );

                word matrixSolverName(fieldName);

                if(fieldTranslationToMatrixSolver_.found(matrixSolverName))
                {
                  matrixSolverName =
                      fieldTranslationToMatrixSolver_[matrixSolverName];
                }

                const dictionary& matrixSolverDict =
                    mesh_.solverDict(matrixSolverName);

                if (debug)
                {
                    Info
                        << "Mx solver "
                        << fieldName
                        << matrixSolverDict
                        << nl << endl;
                }

                for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
                {
                    double normFactor(1.0);

                    word normFactorName = "normFactor";
                    if (cmpt==1)
                        normFactorName += word("Y");
                    else if (cmpt==2)
                        normFactorName += word("Z");
                    else if (cmpt>0)
                        normFactorName += word("UNKNOWN");

                    if (substep == 0)
                    {
                        if (matrixSolverDict.found(normFactorName))
                        {
                            matrixSolverDict.lookup(normFactorName)
                          >> normFactor;
                        }
                    }
                    else if (substep > 0)
                    {
                        const word normFactorNameList(normFactorName + "List");
                        scalarList extraNormFactors;

                        if (matrixSolverDict.found(normFactorNameList))
                        {
                            matrixSolverDict.lookup(normFactorNameList) >>
                                extraNormFactors;

                            if((substep-2) < extraNormFactors.size())
                            {
                                normFactor = extraNormFactors[substep-1];
                            }
                        }
                    }

                    if (debug)
                    {
                        Info<< "Norm " << fieldName
                            << " is " << normFactor
                            << ", ";
                    }

                    if (component(validComponents, cmpt) != -1)
                    {
                        const double value =
                            component(residual, cmpt)*normFactor;

                        if (writeToFile_)
                        {
                            file(0) << tab << value;
                            file(1) << tab << component(residual, cmpt);
                        }

                        Log
                            << fieldName
                              + word(pTraits<Type>::componentNames[cmpt])
                              + subStepName(substep)
                            << token::ASSIGN << token::SPACE
                            << value << token::END_STATEMENT
                            << token::SPACE;
                    }
                }
            }
            else // if(substep >= sp.size())
            {
                typename pTraits<Type>::labelType validComponents
                (
                    mesh_.validComponents<Type>()
                );

                for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
                {
                    if (component(validComponents, cmpt) != -1)
                    {
                        if (writeToFile_)
                        {
                            file(0) << tab << "N/A";
                            file(1) << tab << "N/A";
                        }

                        Log
                            << fieldName
                              + word(pTraits<Type>::componentNames[cmpt])
                              + subStepName(substep)
                            << token::ASSIGN << token::SPACE
                            << "N/A" << token::END_STATEMENT
                            << token::SPACE;
                    }
                }
            }
        }
        else
        {
            if (writeToFile_)
            {
                file(0) << tab << "N/A";
                file(1) << tab << "N/A";
            }

            Log
                << fieldName + subStepName(substep)
                << token::ASSIGN << token::SPACE
                << "N/A" << token::END_STATEMENT
                << token::SPACE;
        }
    }
}


// ************************************************************************* //
