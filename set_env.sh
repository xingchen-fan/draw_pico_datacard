#!/bin/bash

RUN_KERNEL=$(uname -r | cut -d '-' -f1)

# Setting up root using cmssw at ucsb servers. 
# Please setup ROOT yourself if you are not working on a ucsb server.
if [[ $HOSTNAME =~ .*physics.ucsb.edu || $HOSTNAME =~ compute.*.local ]]; then

  . /cvmfs/cms.cern.ch/cmsset_default.sh

  if [ "$RUN_KERNEL" == "3.10.0" ]; then
    export SCRAM_ARCH=slc7_amd64_gcc700
    cd /net/cms29/cms29r0/pico/cc7/CMSSW_10_2_11_patch1/src
  elif [ "$RUN_KERNEL" == "2.6.32" ]; then
    cd /net/cms29/cms29r0/pico/CMSSW_10_2_11_patch1/src
  fi
  
  . /cvmfs/cms.cern.ch/cmsset_default.sh
  eval `scramv1 runtime -sh`
  cd -
fi

if [ -z ${ROOTSYS} ]; then
  echo "Please setup ROOT"
fi

export SCONSFLAGS="-j $(nproc --all)"

source $(dirname $(readlink -e "$BASH_SOURCE"))/modules/jb_utils/set_env.sh
source $(dirname $(readlink -e "$BASH_SOURCE"))/modules/queue_system/set_env.sh

# Example: usernames=("jbkim ana")
usernames=("jbkim oshiro")
if [[ "${usernames[@]}" =~ "$USER" ]]; then
  if [ "$HOSTNAME" = cms37.physics.ucsb.edu ]; then
    export LOCAL_PICO_DIR="/data/"
    echo "LOCAL_PICO_DIR is set to $LOCAL_PICO_DIR. Set it so blank if you want to use net files."
  fi
else
  export LOCAL_PICO_DIR=""
fi
