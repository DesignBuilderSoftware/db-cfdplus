/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's atmBoundaryLayer* incompressible boundary condition
    classes.
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

#include "atmBLVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

atmBLVelocityFvPatchVectorField::
atmBLVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    n_(pTraits<vector>::zero),
    z_(pTraits<vector>::zero),
    Uref_(10),
    Href_(10),
    deltaref_(270),
    aref_(0.14),
    delta_(370),
    a_(0.22),
    zGround_(0.0)
{}


atmBLVelocityFvPatchVectorField::
atmBLVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false),
    n_(dict.lookup("n")),
    z_(dict.lookup("z")),
    Uref_(dict.lookup<scalar>("U_met")),
    Href_(dict.lookup<scalar>("H_met")),
    deltaref_(dict.lookup<scalar>("delta_met")),
    aref_(dict.lookup<scalar>("a_met")),
    delta_(dict.lookup<scalar>("delta")),
    a_(dict.lookup<scalar>("a")),
    zGround_("zGround", dict, p.size())
{
    if (mag(n_) < SMALL)
    {
        FatalErrorInFunction
            << "magnitude of n must be greater than zero"
            << abort(FatalError);
    }

    n_ /= mag(n_);

    const vectorField& c = patch().Cf();
    const scalarField coord(c & z_);
    scalarField Un(coord.size());

    forAll(coord, i)
    {
        Un[i] = Uref_*pow(deltaref_/(Href_+SMALL),aref_)*pow((coord[i]-zGround_[i])/delta_,a_);
    }

    vectorField::operator=(n_*Un);
}


atmBLVelocityFvPatchVectorField::
atmBLVelocityFvPatchVectorField
(
    const atmBLVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    n_(ptf.n_),
    z_(ptf.z_),
    Uref_(ptf.Uref_),
    Href_(ptf.Href_),
    deltaref_(ptf.deltaref_),
    aref_(ptf.aref_),
    delta_(ptf.delta_),
    a_(ptf.a_),
    zGround_(mapper(ptf.zGround_))
{}


atmBLVelocityFvPatchVectorField::
atmBLVelocityFvPatchVectorField
(
    const atmBLVelocityFvPatchVectorField& blpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(blpvf, iF),
    n_(blpvf.n_),
    z_(blpvf.z_),
    Uref_(blpvf.Uref_),
    Href_(blpvf.Href_),
    deltaref_(blpvf.deltaref_),
    aref_(blpvf.aref_),
    delta_(blpvf.delta_),
    a_(blpvf.a_),
    zGround_(blpvf.zGround_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void atmBLVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
    m(zGround_,zGround_);
}


void atmBLVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
    const atmBLVelocityFvPatchVectorField& blptf =
        refCast<const atmBLVelocityFvPatchVectorField>(ptf);

    zGround_.rmap(blptf.zGround_, addr);
}


void atmBLVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "n", n_);
    writeEntry(os, "z", z_);
    writeEntry(os, "U_met", Uref_);
    writeEntry(os, "H_met", Href_);
    writeEntry(os, "delta_met", deltaref_);
    writeEntry(os, "a_met", aref_);
    writeEntry(os, "delta", delta_);
    writeEntry(os, "a", a_);
    writeEntry(os, "zGround", zGround_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    atmBLVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
