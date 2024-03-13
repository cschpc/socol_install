# vim: sw=2
#!/bin/bash

# modules

if [ "1" ]; then
  source artefacts/modules_cray
fi

GERACLIS_ROOT=/project/project_462000470/socol/cray
SOCOL_SRC=$GERACLIS_ROOT/SOCOLv4.0-testing
# YAXT_SRC="..."
# CDI_SRC=
# CDO_SRC=

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
if test -f "$SOCOL_SRC/configure_socol_lumi_cray"; then
  echo "\n********************************************************\n"
  echo "$SOCOL_SRC/configure_socol_lumi_cray exists, overwrite? (y/N)"
  read yesNO
  if [ "${yesNO}" = "y" ]; then
    cp artefacts/configure_socol_lumi_cray $SOCOL_SRC/
  fi
else
  cp artefacts/configure_socol_lumi_cray $SOCOL_SRC/
fi

echo "Copying configure.ac $SOCOL_SRC/src/oasis3mct/configure.ac"
cp artefacts/configure.ac $SOCOL_SRC/src/oasis3mct/configure.ac
pushd $SOCOL_SRC/src/oasis3mct/ && autoconf
popd

# Configure SOCOL
pushd $SOCOL_SRC && GERACLIS_ROOT=$GERACLIS_ROOT ./configure_socol_lumi_cray --oasis --echam --mpiom --prefix=$SOCOL_SRC

# Make oasis
cd $SOCOL_SRC/src/oasis3mct && make -j2 install

# Make SOCOL
cd $SOCOL_SRC && make -j12 install
popd
