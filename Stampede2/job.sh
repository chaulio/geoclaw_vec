#!/bin/bash
#----------------------------------------------------
# Sample SLURM job script
#   for TACC Stampede2 KNL nodes
#
#   *** OpenMP Job on Normal Queue ***
# 
# Last revised: 27 Jun 2017
#
# Notes:
#
#   -- Launch this script by executing
#   -- Copy/edit this script as desired.  Launch by executing
#      "sbatch knl.openmp.slurm" on a Stampede2 login node.
#
#   -- OpenMP codes run on a single node (upper case N = 1).
#        OpenMP ignores the value of lower case n,
#        but slurm needs a plausible value to schedule the job.
#
#   -- Default value of OMP_NUM_THREADS is 1; be sure to change it!
#
#   -- Increase thread count gradually while looking for optimal setting.
#        If there is sufficient memory available, the optimal setting
#        is often 68 (1 thread per core) or 136 (2 threads per core).

#----------------------------------------------------

#SBATCH -J geoclaw           # Job name
#SBATCH -o myjob.o%j       # Name of stdout output file
#SBATCH -e myjob.e%j       # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for OpenMP)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for OpenMP)
#SBATCH -t 01:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=chauliof@gmail.com
#SBATCH --mail-type=end    # Send email at begin and end of job
#---SBATCH -A myproject       # Allocation name (req'd if you have more than 1)
#SBATCH --export num_threads,version,LD_LIBRARY_PATH,CLAW
# Other commands must follow all #SBATCH directives...

LOGGING='tee logs/knl_'$num_threads'_v'$version'.log'

module list | $LOGGING
pwd | $LOGGING
date | $LOGGING

source $HOME/.bashrc

echo CLAW=$CLAW | $LOGGING
echo PATH=$PATH | $LOGGING
echo LD_LIBRARY_PATH=$LD_LIBRARY_PATH | $LOGGING

machine=knl

WORKDIR=$WORK/geoclaw_vec/chile2010
cd $WORKDIR

python setrun.py
mkdir -p logs

export OMP_NUM_THREADS=$num_threads

echo NUM_THREADS=$OMP_NUM_THREADS | $LOGGING

./xgeoclaw-$version.$machine | $LOGGING


