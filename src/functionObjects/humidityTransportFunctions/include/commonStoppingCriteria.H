/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | Developed for blueCFD(R)-Kernel
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2019 FSD blueCAPE Lda
     \\/     M anipulation  |     http://www.bluecape.com.pt
-------------------------------------------------------------------------------
License
    This file was designed to work with OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

bool stoppingCriteria
(
    label iter,
    scalar& resd,
    scalar& resd_0,
    const scalar initialResidual,
    const dictionary& matrixSolverDict,
    const scalar& residualDropTol,
    const scalar& minimumResidual
)
{
    scalar rdrop = 1.0;
    scalar normFactor = 1.0;

    if(matrixSolverDict.found("normFactor"))
    {
        matrixSolverDict.lookup("normFactor") >> normFactor;
    }

    // Initialize the residual monitoring
    if(iter == 0)
    {
        resd = initialResidual*normFactor;
        resd_0 = resd;
    }
    else
    {
        resd = initialResidual*normFactor;
        rdrop = resd/(resd_0 + VSMALL);
    }

    return rdrop <= residualDropTol || resd <= minimumResidual;
}
