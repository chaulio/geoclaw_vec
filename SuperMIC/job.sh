#!/bin/bash
#@ environment = $machine ; $version ; $num_threads
#@ wall_clock_limit = 08:00:00
#@ job_name = geoclaw
#@ job_type = parallel
#@ class = phi
#@ node = 1
#@ node_usage = not_shared
#@ initialdir = $(home)/geoclaw_vec/chile2010
#@ output = logs/log-$(jobid).out
#@ error = logs/log-$(jobid).err
#@ notification=always
#@ notify_user=chauliof@gmail.com
#@ queue


. /etc/profile
. /etc/profile.d/modules.sh


host=`cat $LOADL_HOSTFILE | cut -d- -f1`;
    
WORKDIR=$HOME/geoclaw_vec/chile2010
cd $WORKDIR

python setrun.py
mkdir -p logs

export OMP_NUM_THREADS=$num_threads

device=$host-mic0
if [ "$machine" = "host" ] ; then
	device=$host-ib
fi


LOGGING='tee logs/'$machine'_'$OMP_NUM_THREADS'_v'$version'.log'

pwd | $LOGGING
date | $LOGGING
echo Host = $host-ib | $LOGGING
echo MIC = $host-mic0 | $LOGGING


mpiexec -host $device -n 1 ./xgeoclaw-$version.$machine | $LOGGING
