# vim: sw=2
#!/bin/bash

# configuration
GERACLIS_ROOT=/project/project_462000470
SOCOL_ROOT=$GERACLIS_ROOT/socol/gnu
SOCOL_SRC=$SOCOL_ROOT/SOCOLv4.0-testing
BUILDNAME=testing
CONFNAME=configure_socol_lumi

if [ "1" ]; then
  module use artefacts/modulefiles
  ml socol/$BUILDNAME
fi

# cleanup

# make cleans in SOCOL_SRC

SCRIPTDIR=$(pwd)

pushd $SOCOL_SRC && ls

# squeaky clean
if [ "1" ]; then
  cd $SOCOL_SRC && make clean
  cd $SOCOL_SRC/src/oasis3mct && make clean
  cd $SOCOL_SRC/src/oasis3mct && make clean
  cd $SOCOL_SRC/src/oasis3mct/libpsmile && make clean
  cd $SOCOL_SRC/src/oasis3mct/mct && make clean
fi

popd

# Copy modified configs
if test -f "$SOCOL_SRC/$CONFNAME"; then
  echo "\n********************************************************\n"
  echo "$SOCOL_SRC/$CONFNAME exists, overwrite? (y/N)"
  read yesNO
  if [ "${yesNO}" = "y" ]; then
    cp artefacts/$CONFNAME $SOCOL_SRC/
  fi
else
  cp artefacts/$CONFNAME $SOCOL_SRC/
fi

echo "Copying configure.ac $SOCOL_SRC/src/oasis3mct/configure.ac"
cp artefacts/configure.ac $SOCOL_SRC/src/oasis3mct/configure.ac
pushd $SOCOL_SRC/src/oasis3mct/ && autoconf
popd

# Configure SOCOL
pushd $SOCOL_SRC && SOCOL_ROOT=$SOCOL_ROOT ./$CONFNAME --oasis --echam --mpiom --prefix=$SOCOL_SRC

# Make oasis
cd $SOCOL_SRC/src/oasis3mct && make -j2 install

# Make SOCOL
cd $SOCOL_SRC && make -j12 install
popd
