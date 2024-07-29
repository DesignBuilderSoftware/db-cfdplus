/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Derived from OpenFOAM's scalarTransport function object.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2023 FSD blueCAPE Lda  http://www.bluecape.com.pt
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
    Foam::functionObjects::estimateTemperatureOnHeatFluxPatches

Group
    grpFieldFunctionObjects

Description
    Revise the temperature for patches in the \c Tinput field that use
    boundary conditions that should recalculate their values, e.g. when there
    is heat flux.

    Added \c onConstruction option, because we need it to load right when it
    starts.

    Example of function object specification:
    \verbatim
    reviseTemperatureBCFDK1
    {
        type        estimateTemperatureOnHeatFluxPatches;
        libs        ("libfunctionObjectsBCFDK.so");
        ...

        onConstruction true;
        neverExecute   true;

        Tinput   T;
        Toutput  TempEstimated;

        Tmin     260;
        Tmax     500;
    }
    \endverbatim

Usage
    \table
        Property   | Description             | Required    | Default value
        type       | type name: estimateTemperatureOnHeatFluxPatches | yes  |
        onConstruction | read when construction is complete | yes |
        neverExecute | do not execute at any time | yes |
        Tinput     | Temperature field name to start with       |  no      | T
        Toutput    | Temperature field name for the result      |  yes     |
        Tmin       | If defined, limit the minimum temperature  |  no      |
        Tmax       | If defined, limit the maximum temperature  |  no      |
    \endtable

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    estimateTemperatureOnHeatFluxPatches.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_estimateTemperatureOnHeatFluxPatches_H
#define functionObjects_estimateTemperatureOnHeatFluxPatches_H

#include "fvMeshFunctionObject.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class estimateTemperatureOnHeatFluxPatches Declaration
\*---------------------------------------------------------------------------*/

class estimateTemperatureOnHeatFluxPatches
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Whether to read during construction
        Switch onConstruction_;

        //- Whether to read during execution (which will clean up)
        Switch neverExecute_;

        //- Name of input and output Temperature fields
        word TinputName_;
        word ToutputName_;

        dimensionedScalar Tmin_;
        dimensionedScalar Tmax_;


    // Private Member Functions

        //- Correct the T fields
        void correct();

public:

    //- Runtime type information
    TypeName("estimateTemperatureOnHeatFluxPatches");


    // Constructors

        //- Construct from Time and dictionary
        estimateTemperatureOnHeatFluxPatches
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        estimateTemperatureOnHeatFluxPatches(const estimateTemperatureOnHeatFluxPatches&) = delete;

    //- Destructor
    virtual ~estimateTemperatureOnHeatFluxPatches();


    // Member Functions

        //- Read the estimateTemperatureOnHeatFluxPatches data
        virtual bool read(const dictionary&);

        //- Do not do anything
        virtual bool execute();

        //- Calculate the estimateTemperatureOnHeatFluxPatches
        virtual bool write();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const estimateTemperatureOnHeatFluxPatches&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
