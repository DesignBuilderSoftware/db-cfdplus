/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2015-2016 FSD blueCAPE Lda  http://www.bluecape.com.pt
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

Note
    - Based on OpenFOAM's calculatedFvPatchField. Simply modified to hack
    into the cyclic patch type and it's not meant to be used in fields that
    are part of equations! The objective is to enforce the storage of our
    calculated values.

\*---------------------------------------------------------------------------*/

#include "calculatedDummyCyclicFvPatchField.H"
#include "fvPatchFieldMapper.H"
#include "cyclicFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::calculatedDummyCyclicFvPatchField<Type>::calculatedDummyCyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    calculatedFvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::calculatedDummyCyclicFvPatchField<Type>::calculatedDummyCyclicFvPatchField
(
    const calculatedDummyCyclicFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    calculatedFvPatchField<Type>(ptf, p, iF, mapper)
{
    if (notNull(iF) && mapper.hasUnmapped())
    {
        WarningInFunction
            << "On field " << iF.name() << " patch " << p.name()
            << " patchField " << this->type()
            << " : mapper does not map all values." << nl
            << "    To avoid this warning fully specify the mapping in derived"
            << " patch fields." << endl;
    }
}


template<class Type>
Foam::calculatedDummyCyclicFvPatchField<Type>::calculatedDummyCyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    calculatedFvPatchField<Type>(p, iF, dict, valueRequired)
{}


template<class Type>
Foam::calculatedDummyCyclicFvPatchField<Type>::calculatedDummyCyclicFvPatchField
(
    const calculatedDummyCyclicFvPatchField<Type>& ptf
)
:
    calculatedFvPatchField<Type>(ptf)
{}


template<class Type>
Foam::calculatedDummyCyclicFvPatchField<Type>::calculatedDummyCyclicFvPatchField
(
    const calculatedDummyCyclicFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    calculatedFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type>
void Foam::calculatedDummyCyclicFvPatchField<Type>::write(Ostream& os) const
{
    // We have to override how the original code writes to file, otherwise we
    // risk having patchType being written twice. Which is why we don't call
    // "fvPatchField<Type>::write(os)"

    writeEntry(os, "type", type());

    //This is so that we can enforce the type of underlying patch type
    writeEntry(os, "patchType", cyclicFvPatch::typeName);

    writeEntry(os, "value", *this);
}


// ************************************************************************* //
