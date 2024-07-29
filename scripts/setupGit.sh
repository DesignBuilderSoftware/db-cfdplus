#!/bin/bash

function defaultConfiguration()
{
  git config --local core.autocrlf input
  git config --local core.eol lf
  git config --local core.safecrlf warn

  git config --local color.diff auto
  git config --local color.status auto
  git config --local color.branch auto
  git config --local color.interactive true
}

cd ${0%/*}/.. || exit 1    # run from the parent directory

set -x

defaultConfiguration

if [ -d .git ]
then
  cp scripts/pre-commit .git/hooks/pre-commit
else
  GITPATH=$(sed -e 's=gitdir: ==' .git)
  cp scripts/pre-commit "${GITPATH}"/hooks/pre-commit
fi

set +x

echo "blueCFD-Kernel Code-On-FOAM Git settings deployed."
