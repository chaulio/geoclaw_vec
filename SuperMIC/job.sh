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

export OMP_NUM_THREADS=$num_threads

#env vars for nested OpenMP
#export OMP_NESTED=TRUE
#export OMP_PROC_BIND=spread,close
#export OMP_PLACES=threads
#export OMP_DISPLAY_ENV=verbose
#export KMP_HOT_TEAMS_MODE=1
#export KMP_TEAMS_MAX_LEVEL=2
export KMP_AFFINITY=verbose,$KMP_AFFINITY


RUNDIR="run_"$machine"_"$num_threads"_v"$version
mkdir -p $RUNDIR
cd $RUNDIR
python ../setrun.py


device=$host-mic0
if [ "$machine" = "host" ] ; then
	device=$host-ib
fi


LOGGING='tee '$WORKDIR'/logs/'$machine'_'$OMP_NUM_THREADS'_v'$version'.log'

pwd 2>&1 | $LOGGING
date 2>&1 | $LOGGING
echo Host = $host-ib | $LOGGING
echo MIC = $host-mic0 | $LOGGING


mpiexec -host $device -n 1 ../xgeoclaw-$version.$machine 2>&1 | $LOGGING
