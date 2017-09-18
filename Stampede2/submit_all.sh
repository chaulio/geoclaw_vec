for machine in knl  ; do
	for version in 1 2 ; do
        for num_threads in 1 8 ; do 
	
            export machine
            export version
            export num_threads
            
            echo "==> Submitting: $machine v$version, num_threads: $num_threads"
            sbatch $WORK/geoclaw_vec/Stampede2/job.sh
         done;
	done;
done;
