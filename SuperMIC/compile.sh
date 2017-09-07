WORKDIR=$HOME/geoclaw_vec/chile2010
cd $WORKDIR

flags="-fast -ipo -fopenmp -align array64byte -fp-model=precise "
export FC=ifort

for version in 1 2 3 ; do
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
