#!/bin/bash

# run_benchmarks.sh
# This test script compiles each benchmark problem in serial and parallel
# then executes the parallel program in 2D. The length of the test can be
# controlled by selecting default (1,000 steps), long (5,000 steps), or
# extra long (25,000 steps). The resulting checkpoints are then converted
# to PNG images to spot check problems in the code or environment, including
# execution errors and numerical instabilities.
#
# Questions/comments to trevor.keller@gmail.com (Trevor Keller)

# Valid flags are:
# --noexec  build in serial and parallel, but do not execute
# --noviz   do not convert data for visualization
# --force   pass -B flag to make
# --np X    pass mpirun X ranks (e.g., --np 3 yields mpirun -np 3)
# --short   execute tests .1x default
# --long    execute tests  5x longer  than default
# --extra   execute tests 25x longer  than default
# --clean   delete generated files (binaries, data, imges) after test completes

# Set output colors
RED='\033[0;31m'
GRN='\033[0;32m'
WHT='\033[0m' # No Color

# Initialize timer and completion counters
tstart=$(date +%s)
nSerErr=0
nParErr=0
nRunErr=0
nSerBld=0
nParBld=0
nParRun=0
MFLAG="-s"

# Set execution parameters
ITERS=4000
INTER=250
CORES=4
COREMAX=$(nproc)
if [[ $CORES -gt $COREMAX ]]
then
	CORES=$COREMAX
fi

problems=$(pwd)

# Define the directories to execute, in order. Formatting matters.
exdirs=("periodic/" \
"no-flux/" \
"T-shape/")

exlabels=("a" \
"b" \
"c")

echo -n "Building problems in serial and parallel"

while [[ $# -gt 0 ]]
do
	key="$1"
	case $key in
		--force)
			echo -n ", forcing build"
			MFLAG="-Bs"
		;;
			--clean)
			echo -n ", cleaning up after"
		CLEAN=true
		;;
		--noexec)
			echo -n ", not executing"
			NEXEC=true
		;;
		--short)
			ITERS=$(($ITERS/10))
			INTER=$ITERS
		;;
		--long)
			ITERS=$((10*$ITERS))
			INTER=$(($ITERS/10))
		;;
		--extra)
			ITERS=$((100*$ITERS))
			INTER=$(($ITERS/100))
		;;
		--steadystate)
			ITERS=51200000
			INTER=102400
		;;
		--noviz)
			echo -n ", no PNG output"
			NOVIZ=true
		;;
		--np)
			shift
			CORES=$1
		;;
		*)
			echo "WARNING: Unknown option ${key}."
			echo
    ;;
	esac
	shift # pop first entry from command-line argument list, reduce $# by 1
done

if [[ ! $NEXEC ]]
then
		echo ", taking $ITERS steps, using $CORES/$COREMAX MPI ranks"
else
	echo
fi

if [[ ! $NOVIZ ]]
then
	# Remove existing images first. If the test fails,
	# no images in the directory is an obvious red flag.
	rm -f test.*.png
	if [[ $(which mmsp2png) == "" ]]
	then
		# Consult doc/MMSP.manual.pdf if this fails.
		echo "mmsp2png utility not found. Please check your installation,"
		echo " or pass --noviz flag to suppress PNG output."
	fi
fi

echo "--------------------------------------------------------------------------"

rm -rf ./*/meta.yml ./*/error.log
codeversion=$(git submodule status | awk '{print $1}')
repoversion=$(git rev-parse --verify HEAD)
cpufreq=$(lscpu | grep "max MHz" | awk '{print $NF/1000}')
if [[ $cpufreq == "" ]]
then
	cpufreq=$(grep -m1 MHz /proc/cpuinfo | awk '{print $NF/1000}')
fi
ndots=$(($ITERS / $INTER))
sumspace=$((32 - $ndots/2))

n=${#exdirs[@]}
for (( i=0; i<$n; i++ ))
do
	exstart=$(date +%s)
	j=$(($i+1))
	printf "%s %-${sumspace}s" ${exlabels[$i]} ${exdirs[$i]}
	cd $problems/${exdirs[$i]}

	# Write simulation particulars. Should work on any Debian-flavored GNU/Linux OS.
	echo "---" >>meta.yml
	echo "benchmark_id: 1${exlabels[$i]}" >>meta.yml
	echo "" >>meta.yml
	echo "metadata:" >>meta.yml
	echo "  # Describe the runtime environment" >>meta.yml
	echo "  summary: MPI parallel Travis-CI benchmark with MMSP, ${exdirs[$i]/\//} domain" >>meta.yml
	echo "  author: Trevor Keller" >>meta.yml
	echo "  email: trevor.keller@nist.gov" >>meta.yml
	echo "  date: $(date -R)" >>meta.yml
	echo "  hardware:" >>meta.yml
	echo "    # Required hardware details" >>meta.yml
	echo "    architecture: $(uname -m)" >>meta.yml
	echo "    cores: ${CORES}" >>meta.yml
	echo "    # Optional hardware details" >>meta.yml
	echo "    details:" >>meta.yml
	echo "      - name: clock" >>meta.yml
	echo "        values: ${cpufreq}" >>meta.yml
	echo "        # unit: GHz" >>meta.yml
	echo "  software:" >>meta.yml
	echo "    name: mmsp" >>meta.yml
	echo "    url: https://github.com/mesoscale/mmsp" >>meta.yml
	echo "    version: \"4\"" >>meta.yml
	echo "    repo:" >>meta.yml
	echo "      url: https://github.com/mesoscale/mmsp/tree/develop" >>meta.yml
	echo "      version: \"${codeversion}\"" >>meta.yml
	echo "  implementation:" >>meta.yml
	echo "    end_condition: time limit, Travis CI runs die after 50 minutes total" >>meta.yml
	echo "    repo:" >>meta.yml
	echo "      url: https://github.com/mesoscale/MMSP-spinodal-decomposition-benchmark/tree/master/${exdirs[$i]}" >>meta.yml
	echo "      version: \"${repoversion}\"" >>meta.yml
	echo "      # badge: https://api.travis-ci.org/mesoscale/MMSP-spinodal-decomposition-benchmark.svg?branch=master" >>meta.yml
	echo "    details:" >>meta.yml
	echo "      - name: mesh" >>meta.yml
	echo "        values: uniform rectilinear" >>meta.yml
	echo "      - name: numerical_method" >>meta.yml
	echo "        values: explicit finite difference" >>meta.yml
	echo "      - name: compiler" >>meta.yml
	echo "        values: GNU mpic++" >>meta.yml
	echo "      - name: parallel_model" >>meta.yml
	echo "        values: MPI" >>meta.yml
	echo "" >>meta.yml

	if make $MFLAG
	then
		((nSerBld++))
	else
		((nSerErr++))
	fi
	if make $MFLAG parallel
	then
		((nParBld++))
	else
		((nParErr++))
	fi
	if [[ -f parallel ]] && [[ ! $NEXEC ]]
	then
		# Run the example in parallel, for speed.
		rm -f test.*.dat
		(/usr/bin/time -f "  - name: run_time\n    values: {\"time\": %e, \"unit\": seconds}\n  - name: memory_usage\n    values: {\"value\": %M, \"unit\": KB}" bash -c \
		"/usr/bin/mpirun.openmpi -np $CORES ./parallel --example 2 test.0000.dat 1>>meta.yml 2>>error.log && \
		/usr/bin/mpirun.openmpi -np $CORES ./parallel test.0000.dat $ITERS $INTER 1>>meta.yml 2>>error.log") &>>meta.yml &

		# Travis CI quits after 10 minutes with no CLI activity. Give it an indication that things are running.
		JOBID=$!
		sleep 30
		OLDFILES=$(ls -1 test*.dat | wc -l)
		while kill -0 "$JOBID" &>/dev/null
		do
			sleep 14
			NEWFILES=$(ls -1 test*.dat | wc -l)
			if [[ $NEWFILES > $OLDFILES ]]
			then
				# A checkpoint was written while we slept. Tell the terminal.
				echo -n '•'
				OLDFILES=$NEWFILES
			fi
		done

		# Return codes are not reliable. Save errors to disk for postmortem.
		if [[ -f error.log ]] && [[ $(wc -w error.log) > 1 ]]
		then
			for i in `seq 1 ${sumspace}`;
			do
				echo -n " "
			done
			echo -e "${RED} --FAILED--${WHT}"
			((nRunErr++))
			if [[ -f error.log ]]
			then
				echo "      error.log has the details (head follows)"
				head error.log | sed -e 's/^/      /'
			fi
		else
			((nParRun++))
			rm -f error.log
			if [[ ! $NOVIZ ]]
			then
				# Show the result
				for f in *.dat
				do
					mmsp2png --zoom $f 1>/dev/null 2>>error.log
				done
			fi
			exfin=$(date +%s)
			exlapse=$(echo "$exfin-$exstart" | bc -l)
			printf "${GRN}%${sumspace}d seconds${WHT}\n" $exlapse
		fi
	else
		exfin=$(date +%s)
		exlapse=$(echo "$exfin-$exstart" | bc -l)
		printf "${GRN}%${sumspace}d seconds${WHT}\n" $exlapse
	fi
	# Clean up binaries and images
	if [[ $CLEAN ]]
	then
		make -s clean
		rm -f test.*.dat
		rm -f error.log
		rm -f test.*.png
	fi
done

cd ${problems}

tfinish=$(date +%s)
elapsed=$(echo "$tfinish-$tstart" | bc -l)
echo "--------------------------------------------------------------------------"
printf "Elapsed time: %52d seconds\n" $elapsed
echo
printf "%2d serial   problems compiled successfully" $nSerBld
if [[ $nSerErr > 0 ]]
then
	printf ", %2d failed" $nSerErr
fi
printf "\n%2d parallel problems compiled successfully" $nParBld
if [[ $nParErr > 0 ]]
then
	printf ", %2d failed" $nParErr
fi
printf "\n%2d parallel problems executed successfully" $nParRun
if [[ $nRunErr > 0 ]]
then
	printf ", %2d failed" $nRunErr
fi
echo
cd ${problems}

AllERR=$(echo "$nSerErr+$nParErr+$nRunErr" | bc -l)
if [[ $AllERR > 0 ]]
then
	echo
	echo "Build error(s) detected: ${AllERR} tests failed."
	exit 1
fi
