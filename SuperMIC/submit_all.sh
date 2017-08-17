for machine in mic host ; do
	for version in 1 2 ; do
        for openmp in yes no ; do 
	
            export machine
            export version
            export openmp
            
            echo "==> Submitting: $machine $version, openmp: $openmp"
            llsubmit /home/hpc/pr45fi/di49rew/geoclaw_vec/SuperMIC/execute.sh
         done;
	done;
done;