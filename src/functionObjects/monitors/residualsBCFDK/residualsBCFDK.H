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

Class
    Foam::functionObjects::residualsBCFDK

Description
    Writes out the initial residual for specified fields.

    Example of function object specification:
    \verbatim
    residualsBCFDK
    {
        type            residualsBCFDK;
        writeControl   timeStep;
        writeInterval  1;
        debug          0;
        writeToFile    true; //false to override and don't write to file
        log            true; //show on screen
        fields
        (
            U
            p
            p
            Contaminant_MassRatio
        );
        substeps  //optional
        (
            0
            0
            1
        );
        fieldToSolver //optional
        (
            Contaminant_MassRatio Contaminant
        );
    }
    \endverbatim

    Output data is written to the dir postProcessing/residualsBCFDK/\<timeDir\>/
    For vector/tensor fields, e.g. U, where an equation is solved for each
    component, the largest residual of each component is written out.
    If a field is solved more than once per time step, then explicitly state
    which sub-step you want to write, where 0 is the first one (default).

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles
    Foam::functionObjects::timeControl

SourceFiles
    residualsBCFDK.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_residualsBCFDK_H
#define functionObjects_residualsBCFDK_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class residualsBCFDK Declaration
\*---------------------------------------------------------------------------*/

class residualsBCFDK
:
    public fvMeshFunctionObject,
    public logFiles
{
protected:

    // Protected data

        //- Fields to write residuals
        wordList fieldSet_;

        //- Respective substeps for each field
        labelList fieldSetSubStep_;

        //- Corresponding table of field names to matrix solver names
        HashTable<word> fieldTranslationToMatrixSolver_;

        //- should it really write to file
        bool writeToFile_;


    // Protected Member Functions

        //- Output field header information
        template<class Type>
        void writeFileHeader
        (
            const label i,
            const word& fieldName,
            const label substep
        );

        //- Output file header information
        virtual void writeFileHeader(const label i);

        //- Calculate the field min/max
        template<class Type>
        void writeResidual(const word& fieldName, const label substep);

        word subStepName(const label substep);

public:

    //- Runtime type information
    TypeName("residualsBCFDK");


    // Constructors

        //- Construct from Time and dictionary
        residualsBCFDK
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        residualsBCFDK(const residualsBCFDK&) = delete;


    //- Destructor
    virtual ~residualsBCFDK();


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the residuals
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const residualsBCFDK&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "residualsBCFDKTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
