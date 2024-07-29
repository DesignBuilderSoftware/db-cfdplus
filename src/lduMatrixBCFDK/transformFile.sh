#!/bin/bash

INPUTFILE="${1}"
KEYWORD=""

for keyword in GAMGSolver PBiCGStab
do
    if [ "${INPUTFILE/$keyword/}" != "${INPUTFILE}" ]
    then
        KEYWORD=$keyword
    fi
done

if [ -n "${KEYWORD}" ]
then
    sed -e "s=${KEYWORD}=BCFDK_${KEYWORD}=g" \
        -e '/scalar normFactor/a\
    \
    \/\/FSD blueCAPE Lda 2018\/2023: store the normFactor\
    #include "normFactorStore.H"' \
        -e 's/"GAMG"/"BCFDK_GAMG"/' \
        -e '/public lduMatrix::solver/{N;s/$/\n\
    \/\/FSD blueCAPE Lda 2018\/2023: Keep a copy of the solver controls\
    const dictionary\& solverControls_;\n/}' \
        -e '/^        solverControls/{N;s/\()\)\(,*\)$/\1,\
    solverControls_(solverControls)\2/}' \
        "${INPUTFILE}"
fi

# ----------------------------------------------------------------- end-of-file
