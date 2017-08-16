#!/bin/bash
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

for machine in mic host ; do
	for version in 1 2 ; do
		if [ "$machine" = "host" ]
		then
            device=$host-ib
			#./xgeoclaw-$version.$machine > logs/log-$version.$machine
		else
            device=$host-mic0
		fi
		
		mpiexec -host $device -n 1 ./xgeoclaw-$version.$machine | tee logs/log-$version.$machine
		
	done;
done;