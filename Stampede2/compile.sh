WORKDIR=$WORK/geoclaw_vec/chile2010
cd $WORKDIR

echo $CLAW

export machine=knl

export FFLAGS="-O3 -ipo -fopenmp -xMIC-AVX512 -align array64byte -fp-model=precise "
export FC=ifort

for version in 1 2  ; do
		export VERSION=$version
		echo "==> Compiling version $version for $machine with flags: $FFLAGS"
		make new
		mv xgeoclaw xgeoclaw-$version.$machine
		echo "DONE: xgeoclaw-$version.$machine"
done;
