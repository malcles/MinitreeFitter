#!/bin/bash

system=`uname -r`

if [ ! -e Makefile ]; then
    . MTools/setupmtools.sh
fi

if    echo $system | grep -q "el6" ; then
    echo "==> SLC6 detected"
    slcCompilerOpt=x86_64-slc6-gcc47-opt
elif  echo $system | grep -q "el5" ; then
    echo "==> SLC5 detected"
    slcCompilerOpt=x86_64-slc5-gcc47-opt
else
    echo "==> unknown operating system"
fi

#source etc/scripts/env.sh
ROOTSYS=/afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.07a/${slcCompilerOpt}/root/
PYTHONDIR=/afs/cern.ch/sw/lcg/external/Python/2.7.3/${slcCompilerOpt}/
XRDCP=/afs/cern.ch/sw/lcg/external/xrootd/3.2.7/${slcCompilerOpt}/
VDTDIR=/afs/cern.ch/sw/lcg/external/vdt/0.3.2/${slcCompilerOpt}/
. /afs/cern.ch/sw/lcg/external/gcc/4.7.2/${slcCompilerOpt}/setup.sh
. $ROOTSYS/bin/thisroot.sh

xAnaExtLIBS=etc/external/
export VDTDIR

export PYTHONPATH=$PYTHONDIR:$ROOTSYS/lib
LD_LIBRARY_PATH=$ROOTSYS/lib:$PYTHONDIR/lib:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=`pwd`/lib:$XRDCP/lib64:$LD_LIBRARY_PATH:.
LD_LIBRARY_PATH=$xAnaExtLIBS/lib/:$LD_LIBRARY_PATH

export LD_LIBRARY_PATH
export PATH=`pwd`/bin/:`pwd`/etc/scripts/:$PYTHONPATH/bin/:$XRDCP/bin/:$PATH


