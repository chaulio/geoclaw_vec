WORKDIR=$HOME/geoclaw_vec/chile2010
cd $WORKDIR

# if these vars are not set, give them these default values:
: ${FC:=ifort}
: ${flags:="-fast -ipo -fopenmp -align array64byte -fp-model=precise "}
: ${versions:="0 1 2"}
: ${machines:="mic host"}


for version in $versions ; do
	for machine in $machines ; do
		if [ "$machine" = "mic" ]
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
