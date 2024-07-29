/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2014-2018 FSD blueCAPE Lda  http://www.bluecape.com.pt
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

#include "flowRateInletDirectionVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "one.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class RhoType>
void Foam::flowRateInletDirectionVelocityFvPatchVectorField::updateValues
(
    const RhoType& rho
)
{
    const scalar t = db().time().timeOutputValue();

    if (extrapolateProfile_)
    {
        //WARNING: This part of the code was modified and is not yet tested!!

        vectorField Up(this->patchInternalField());

        // Patch normal extrapolated velocity
        scalarField nUp(n_ & Up);

        // Remove the normal component of the extrapolate patch velocity
        Up -= nUp*n_;

        // Remove any reverse flow
        nUp = min(nUp, scalar(0));

        const scalar flowRate = flowRate_->value(t);
        const scalar estimatedFlowRate = -gSum(rho*(this->patch().magSf()*nUp));

        if (estimatedFlowRate/flowRate > 0.5)
        {
            nUp *= (mag(flowRate)/mag(estimatedFlowRate));
        }
        else
        {
            nUp -= ((flowRate - estimatedFlowRate)/gSum(rho*patch().magSf()));
        }

        // Add the corrected normal component of velocity to the patch velocity
        Up += nUp*n_;

        // Correct the patch velocity
        this->operator==(Up);
    }
    else
    {
        const scalar avgU = flowRate_->value(t)/gSum(rho*patch().magSf());
        operator==(avgU*nq_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowRateInletDirectionVelocityFvPatchVectorField::
flowRateInletDirectionVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    n_(p.size()),
    nq_(p.size()),
    flowRate_(),
    volumetric_(false),
    rhoName_("rho"),
    rhoInlet_(0.0),
    extrapolateProfile_(false)
{}


Foam::flowRateInletDirectionVelocityFvPatchVectorField::
flowRateInletDirectionVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, false),
    n_("n", dict, p.size()),
    rhoInlet_(dict.lookupOrDefault<scalar>("rhoInlet", -VGREAT)),
    extrapolateProfile_
    (
        dict.lookupOrDefault<Switch>("extrapolateProfile", false)
    )
{
    if (dict.found("volumetricFlowRate"))
    {
        volumetric_ = true;
        flowRate_ = Function1<scalar>::New("volumetricFlowRate", dict);
        rhoName_ = "rho";
    }
    else if (dict.found("massFlowRate"))
    {
        volumetric_ = false;
        flowRate_ = Function1<scalar>::New("massFlowRate", dict);
        rhoName_ = word(dict.lookupOrDefault<word>("rho", "rho"));
    }
    else
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Please supply either 'volumetricFlowRate' or"
            << " 'massFlowRate' and 'rho'" << exit(FatalIOError);
    }

    //FIXME: nq_ should be computed whenever the mesh is updated!
    //Normalize the direction vector
    nq_ = n_ / mag(n_);

    //Scale the direction vector to be on the same scale on the direction
    //of the patch normal
    nq_ = nq_ / (nq_ & -patch().nf());

    // Value field require if mass based
    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        evaluate(Pstream::commsTypes::blocking);
    }
}


Foam::flowRateInletDirectionVelocityFvPatchVectorField::
flowRateInletDirectionVelocityFvPatchVectorField
(
    const flowRateInletDirectionVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    n_(mapper(ptf.n_)),
    nq_(mapper(ptf.nq_)),
    flowRate_(ptf.flowRate_, false),
    volumetric_(ptf.volumetric_),
    rhoName_(ptf.rhoName_),
    rhoInlet_(ptf.rhoInlet_),
    extrapolateProfile_(ptf.extrapolateProfile_)
{}


Foam::flowRateInletDirectionVelocityFvPatchVectorField::
flowRateInletDirectionVelocityFvPatchVectorField
(
    const flowRateInletDirectionVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    n_(ptf.n_),
    nq_(ptf.nq_),
    flowRate_(ptf.flowRate_, false),
    volumetric_(ptf.volumetric_),
    rhoName_(ptf.rhoName_),
    rhoInlet_(ptf.rhoInlet_),
    extrapolateProfile_(ptf.extrapolateProfile_)
{}


Foam::flowRateInletDirectionVelocityFvPatchVectorField::
flowRateInletDirectionVelocityFvPatchVectorField
(
    const flowRateInletDirectionVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    n_(ptf.n_),
    nq_(ptf.nq_),
    flowRate_(ptf.flowRate_, false),
    volumetric_(ptf.volumetric_),
    rhoName_(ptf.rhoName_),
    rhoInlet_(ptf.rhoInlet_),
    extrapolateProfile_(ptf.extrapolateProfile_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::flowRateInletDirectionVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField::autoMap(m);

    m(n_,n_);
    m(nq_,nq_);
}


void Foam::flowRateInletDirectionVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField::rmap(ptf, addr);

    const flowRateInletDirectionVelocityFvPatchVectorField& tiptf =
        refCast<const flowRateInletDirectionVelocityFvPatchVectorField>
        (ptf);

    n_.rmap(tiptf.n_, addr);
    nq_.rmap(tiptf.nq_, addr);
}


void Foam::flowRateInletDirectionVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (volumetric_ || rhoName_ == "none")
    {
        updateValues(one());
    }
    else
    {
        // Mass flow-rate
        if (db().foundObject<volScalarField>(rhoName_))
        {
            const fvPatchField<scalar>& rhop =
                patch().lookupPatchField<volScalarField, scalar>(rhoName_);

            updateValues(rhop);
        }
        else
        {
            // Use constant density
            if (rhoInlet_ < 0)
            {
                FatalErrorInFunction
                    << "Did not find registered density field " << rhoName_
                    << " and no constant density 'rhoInlet' specified"
                    << exit(FatalError);
            }

            updateValues(rhoInlet_);
        }
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::flowRateInletDirectionVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    writeEntry(os, "n", n_) ;
    writeEntry(os, flowRate_());
    if (!volumetric_)
    {
        writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
        writeEntryIfDifferent<scalar>(os, "rhoInlet", -VGREAT, rhoInlet_);
    }
    writeEntry(os, "extrapolateProfile", extrapolateProfile_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       flowRateInletDirectionVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
