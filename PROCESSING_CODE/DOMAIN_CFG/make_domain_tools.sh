#!/bin/bash

mv $TDIR/DOMAINcfg/src/domain.f90 $TDIR/DOMAINcfg/src/domain.f90.original_off
mv $TDIR/DOMAINcfg/src/dom_oce.f90 $TDIR/DOMAINcfg/src/dom_oce.f90.original_off
mv $TDIR/DOMAINcfg/src/domzgr.f90 $TDIR/DOMAINcfg/src/domzgr.f90.original_off

cp $REPO/SRC/DOMAIN/domain.f90 $TDIR/DOMAINcfg/src/domain.f90
cp $REPO/SRC/DOMAIN/dom_oce.f90 $TDIR/DOMAINcfg/src/dom_oce.f90
cp $REPO/SRC/DOMAIN/domzgr.f90 $TDIR/DOMAINcfg/src/domzgr.f90

module swap craype-network-ofi craype-network-ucx
module swap cray-mpich cray-mpich-ucx
module load cray-hdf5-parallel/1.12.2.1
module load cray-netcdf-hdf5parallel/4.9.0.1

cd $TDIR
./maketools -m archer2_cray -n DOMAINcfg


