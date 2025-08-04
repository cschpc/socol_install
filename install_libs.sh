#!/usr/bin/env bash
source settings.sh

echo SOCOL_SRC=$SOCOL_SRC
source artefacts/modules_gnu

pushd downloads/libcdi-cdi-1.7.2

CXX=CC FC=ftn CC=cc ./configure --prefix=$GERACLIS_ROOT/usr --with-netcdf=$NETCDF_DIR
make -j 12 install

popd

pushd downloads/cdo-2.2.0

CXX=CC FC=ftn CC=cc ./configure --prefix=$GERACLIS_ROOT/usr

make -j 12 install

popd


pushd downloads/yaxt-0.9.1

CXX=CC FC=ftn CC=cc ./configure --prefix=$GERACLIS_ROOT/usr --without-regard-for-quality
make -j 12 install

popd
