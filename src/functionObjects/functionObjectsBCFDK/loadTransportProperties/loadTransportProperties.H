/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Is based on OpenFOAM's readFields function object.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2020 FSD blueCAPE Lda  http://www.bluecape.com.pt
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
    Foam::functionObjects::loadTransportProperties

Description
    Loads the \c transportProperties file onto the registry. Useful when we
    still need incompressible data on compressible flow cases.

    Example of function object specification:
    \verbatim
    loadTransportPropertiesBCFDK
    {
        type        loadTransportProperties;
        libs        ("libfunctionObjectsBCFDK.so");
    }

    See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    readFields.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_loadTransportProperties_H
#define functionObjects_loadTransportProperties_H

#include "fvMeshFunctionObject.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class loadTransportProperties Declaration
\*---------------------------------------------------------------------------*/

class loadTransportProperties
:
    public fvMeshFunctionObject
{
protected:

    // Protected data

        // The transportProperties dictionary
        IOdictionary transportPropertiesDict_;

public:

    //- Runtime type information
    TypeName("loadTransportProperties");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        loadTransportProperties
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        loadTransportProperties(const loadTransportProperties&) = delete;


    //- Destructor
    virtual ~loadTransportProperties();

    // Member Functions

        //- Do nothing
        virtual bool execute();

        //- Do nothing
        virtual bool write();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const loadTransportProperties&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
