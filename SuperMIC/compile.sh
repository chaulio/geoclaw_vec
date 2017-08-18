WORKDIR=$HOME/geoclaw_vec/chile2010
cd $WORKDIR

flags="-O2 -fopenmp -align array64byte -fp-model=precise "
export FC=mpiifort

for version in 1 2 ; do
	for mic in yes no ; do
		if [ "$mic" = "yes" ]
		then
			export FFLAGS="$flags -mmic"
			machine=mic
		else
			export FFLAGS=$flags
			machine=host			
		fi
		export VERSION=$version
		echo "==> Compiling version $version for $machine with flags: $FFLAGS"
		make new
		mv xgeoclaw xgeoclaw-$version.$machine
		echo "DONE: xgeoclaw-$version.$machine"
	done;
done;
