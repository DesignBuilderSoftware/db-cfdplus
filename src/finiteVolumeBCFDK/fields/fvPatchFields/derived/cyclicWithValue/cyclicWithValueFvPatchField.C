/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2015-2018 FSD blueCAPE Lda  http://www.bluecape.com.pt
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

#include "cyclicWithValueFvPatchField.H"
#include "fvPatchFieldMapper.H"
#include "cyclicFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::cyclicWithValueFvPatchField<Type>::cyclicWithValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicFvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::cyclicWithValueFvPatchField<Type>::cyclicWithValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    cyclicFvPatchField<Type>(p, iF, dict)
{
    fvPatchField<Type>::operator=
    (
        Field<Type>("value", dict, p.size())
    );
}


template<class Type>
Foam::cyclicWithValueFvPatchField<Type>::cyclicWithValueFvPatchField
(
    const cyclicWithValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    cyclicFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::cyclicWithValueFvPatchField<Type>::cyclicWithValueFvPatchField
(
    const cyclicWithValueFvPatchField<Type>& ptf
)
:
    cyclicFvPatchField<Type>(ptf)
{}


template<class Type>
Foam::cyclicWithValueFvPatchField<Type>::cyclicWithValueFvPatchField
(
    const cyclicWithValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type>
void Foam::cyclicWithValueFvPatchField<Type>::write(Ostream& os) const
{
    // We have to override how the original code writes to file, otherwise we
    // risk having patchType being written twice. Which is why we don't call
    // "fvPatchField<Type>::write(os)"

    writeEntry(os, "type", type());

    //Essentially it's a direct hack
    writeEntry(os, "patchType", cyclicFvPatch::typeName);

    writeEntry(os, "value", *this);
}


// ************************************************************************* //
