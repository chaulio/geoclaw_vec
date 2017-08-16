#!/bin/bash
#@ environment = $machine ; $version ; $openmp
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

pwd

host=`cat $LOADL_HOSTFILE | cut -d- -f1`;
    
echo Host = $host-ib
echo MIC = $host-mic0

WORKDIR=$HOME/geoclaw_vec/chile2010
cd $WORKDIR

python setrun.py
mkdir -p logs

ssh $host-mic0
pwd
cd $WORKDIR
pwd


if [ "$machine" = "host" ]
then
    device=$host-ib
    export OMP_NUM_THREADS=16
else
    device=$host-mic0
    export OMP_NUM_THREADS=240
fi

if [ "$openmp" = "no" ]
then
    threading=singlethread
    export OMP_NUM_THREADS=1
else
    threading=multithread
fi

mpiexec -host $device -n 1 ./xgeoclaw-$version.$machine | tee logs/log-$version-$threading.$machine
