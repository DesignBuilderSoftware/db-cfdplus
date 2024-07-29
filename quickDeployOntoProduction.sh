#!/bin/sh

cd ${0%/*} || exit 1    # run from this directory

if [ $WM_ARCH_OPTION -eq 64 ]; then

  export TARGET_USER_PLATFORM=../../production/x64

elif [ $WM_ARCH_OPTION -eq 32 ]; then

  export TARGET_USER_PLATFORM=../../production/Win32

fi

OFVERSION=$(wmakePrintBuild -major)
BLUECFD_CORE_NAME="blueCFD-Core"

export TARGET_PLATFORM_FOLDER=$TARGET_USER_PLATFORM/$BLUECFD_CORE_NAME/ofuser-of$OFVERSION/platforms/$WM_OPTIONS

cp -uvr $FOAM_USER_APPBIN/* $TARGET_PLATFORM_FOLDER/bin/
cp -uvr $FOAM_USER_LIBBIN/* $TARGET_PLATFORM_FOLDER/lib/

# ----------------------------------------------------------------- end-of-file
