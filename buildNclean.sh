#!/bin/sh

cd ${0%/*} || exit 1    # run from this directory

set -x
echo "----------------- First Pass -----------------"
./Allwmake

echo "----------------- Second Pass -----------------"
./Allwmake

echo "----------------- CleanUp Pass -----------------"
./Allwclean

# ----------------------------------------------------------------- end-of-file
