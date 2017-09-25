
for machine in host mic
do
	if [ "$machine" = "mic" ] ; then
		num_threads_list="001 002 004 008 015 030 060 120"
	else
		num_threads_list="01 02 04 08 16 32"
	fi

	for version in 1 2
	do
		for num_threads in $num_threads_list
		do
			export machine
			export version
			export num_threads
		        echo "==> Submitting: $machine v$version, threads: $num_threads"
		        llsubmit /home/hpc/pr45fi/di49rew/geoclaw_vec/SuperMIC/job.sh
		done;
	done;
done;



