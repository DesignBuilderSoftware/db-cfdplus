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
    Copyright (C) 2013-2019 FSD blueCAPE Lda  http://www.bluecape.com.pt
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
    Foam::functionObjects::scalarTransportBCFDK

Group
    grpSolversFunctionObjects

Description
    Evolves a passive scalar transport equation.

    - To specify the field name set the \c field entry
    - To employ the same numerical schemes as another field set
      the \c schemesField entry,
    - A constant diffusivity may be specified with the \c D entry,

    - Alternatively if a turbulence model is available a turbulent diffusivity
      may be constructed from the laminar and turbulent viscosities using the
      optional diffusivity coefficients \c alphaD and \c alphaDt (which default
      to 1):
      \verbatim
          D = alphaD*nu + alphaDt*nut
      \endverbatim

    - Beyond what it's done in OpenFOAM by default:

      - Stopping criteria includes \c residualDropTol and \c minimumResidual,
        along with using the normalization factor from the matrix solver, so
        that we can ensure whether the residuals are really dropping or not.

      - A minimum value can be defined, for not allowing values below that
        minimum.

      - The equation can be reset with a specific value, for restarting the
        solution of the equation.

      - Changed the relaxation from the equation to the field.

      - D can be specified as a field that can be used from existing
        on-memory object.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    scalarTransportBCFDK.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_scalarTransportBCFDK_H
#define functionObjects_scalarTransportBCFDK_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "fvOptionList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class scalarTransportBCFDK Declaration
\*---------------------------------------------------------------------------*/

class scalarTransportBCFDK
:
    public fvMeshFunctionObject
{
public:

    // Public data

        //- Diffusion types
        enum diffusionTypes
        {
            //a constant, uniform D_ is specified
            constant,

            // indicate to use existing diffusion coefficient field, which
            // requires that the field name is "Dscalar", where scalar is
            // the name of the scalar field
            existing,

            // automatically inferred from the turbulence model and when that
            // is not available, initialize it as a zero field
            inferred
        };

        //- Ventilator flow directions names
        static const NamedEnum<diffusionTypes, 3> diffusionTypeNames_;

private:

    // Private Data

        //- Name of field to process
        word fieldName_;

        //- Name of flux field (optional)
        word phiName_;

        //- Name of density field (optional)
        word rhoName_;

        //- Diffusion coefficient (optional)
        scalar D_;

        //- Diffusion type to be used
        diffusionTypes dType_;

        //- Minimum allowed value (optional)
        scalar minFieldValue_;

        //- Flag to indicate whether "minValue" was defined
        bool minFieldValueIsSet_;

        //- Reset value before starting to solve (optional)
        scalar resetFieldValue_;

        //- Flag to indicate whether "resetValue" was defined
        bool resetFieldValueIsSet_;

        //- Laminar diffusion coefficient (optional)
        scalar alphaD_;

        //- Turbulent diffusion coefficient (optional)
        scalar alphaDt_;

        //- Number of corrector iterations (optional)
        label nCorr_;

        //- Residual drop tolerance (optional)
        scalar residualDropTol_;

        //- Minimum residual value (optional)
        scalar minimumResidual_;

        //- Name of field whose schemes are used (optional)
        word schemesField_;

        //- Run-time selectable finite volume options, e.g. sources, constraints
        fv::optionList fvOptions_;

        //- The scalar field origin:
        //  - true => create new
        //  - false => look for existing one with the desired scalar field name
        bool createNewScalar_;

    // Private Member Functions

        //- Return the diffusivity field
        tmp<volScalarField> D(const surfaceScalarField& phi) const;

        //- Check the additional stopping criteria
        //  returns true if it should stop.
        bool stoppingCriteria
        (
            label iter,
            scalar& resd,
            scalar& resd_0,
            const scalar initialResidual
        );


public:

    //- Runtime type information
    TypeName("scalarTransportBCFDK");


    // Constructors

        //- Construct from Time and dictionary
        scalarTransportBCFDK
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        scalarTransportBCFDK(const scalarTransportBCFDK&) = delete;

    //- Destructor
    virtual ~scalarTransportBCFDK();


    // Member Functions

        //- Read the scalarTransportBCFDK data
        virtual bool read(const dictionary&);

        //- Calculate the scalarTransportBCFDK
        virtual bool execute();

        //- Do nothing.
        //  The volScalarField is registered and written automatically
        virtual bool write();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const scalarTransportBCFDK&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
