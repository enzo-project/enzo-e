#!/bin/bash

# ==============================
# Verify input parameters
# ==============================

if [[ ("x$1" == "x") || ("x$3" != "x") ]]; then
    echo "Usage: $(basename $0) <Enzo-E output file> [ output directory ]"
    exit 1
fi

# ==============================
# Initialize path and file variables
#==============================

# Get absolute path to performance directory
file=$0
while [[ -L $file ]]; do
      file=`readlink $file`
done
topdir=$(dirname $file)
topdir="$PWD/$topdir"

bindir="`dirname $0`"
if [ "$bindir" == "${bindir#/}" ]; then
    bindir="$PWD/$bindir"
fi
echo "bindir = $bindir"

# Get input file $input
input="$PWD/$1"

# Get output directory $outdir (EZPerf by default)
outdir="$2"
if [ "x$outdir" == "x" ]; then
    outdir="EZPerf"
fi

# ==============================
# Create the output directory
# ==============================

if [[ -e "$outdir" ]]; then
    echo "Warning deleting previous $outdir!"
    rm -rf $outdir
fi
mkdir $outdir
cd $outdir

ADAPT=`awk '/perf:region adapt/{print $(NF-2)}' $input | sort | uniq`
REFRESH=`awk '/perf:region refresh/{print $(NF-2)}' $input | sort | uniq`
METHOD=`awk '/perf:region method/{print $(NF-2)}' $input | sort | uniq`
SOLVER=`awk '/perf:region solver/{print $(NF-2)}' $input | sort | uniq`
MEMORY=`awk '/perf:region cycle /{print $(NF-1)}' $input | sort | uniq`
BALANCE_MAX=`awk '/perf:balance max-/{print $(NF-1)}' $input | sort | uniq`
BALANCE_EFF=`awk '/perf:balance eff-/{print $(NF-2)}' $input | sort | uniq`
MESH=`awk '/perf:mesh /{print $(NF-1)}' $input | sort | uniq`

num_procs=`awk '/CkNumPes/  {print $5}' $input`
num_nodes=`awk '/CkNumNodes/{print $5}' $input`

input_clean="input-clean.data"
grep -v WARNING $input > $input_clean
# ==============================
# Generate data files
# ==============================

echo "nodes $num_nodes procs $num_procs"
if [[ ! -e "cycle.data" ]]; then
    awk '/Simulation cycle/{if ($NF==0) {t0=$2}; print $NF,'"$num_procs"'*($2-t0)}' $input_clean > cycle.data
fi


for adapt in $ADAPT; do
    if [[ ! -e "$adapt.data" ]]; then
           echo "Generating $adapt.data"
           awk '/Simulation cycle/{c=$NF}; /perf:region '"$adapt"' /{print c,$NF}' $input_clean > $adapt.data
       fi
done

for refresh in $REFRESH; do
    if [[ ! -e "$refresh.data" ]]; then
           echo "Generating $refresh.data"
           awk '/Simulation cycle/{c=$NF}; /perf:region '"$refresh"' /{print c,$NF}' $input_clean > $refresh.data
       fi
done

for method in $METHOD; do
    if [[ ! -e "$method.data" ]]; then
        echo "Generating $method.data"
        awk '/Simulation cycle/{c=$NF}; /perf:region '"$method"' /{print c,$NF}' $input_clean > $method.data
    fi
done

for solver in $SOLVER; do
    if [[ ! -e "$solver.data" ]]; then
        echo "Generating $solver.data"
        awk '/Simulation cycle/{c=$NF}; /perf:region '"$solver"' /{print c,$NF}' $input_clean > $solver.data
    fi
done

for memory in $MEMORY; do
    if [[ ! -e "$memory.data" ]]; then
        echo "Generating $memory.data"
        awk '/Simulation cycle/{c=$NF}; /perf:region cycle '"$memory"' /{print c,$NF}' $input_clean > $memory.data
    fi
done
for mesh in $MESH; do
    if [[ ! -e "$mesh.data" ]]; then
        echo "Generating $mesh.data"
        awk '/Simulation cycle/{c=$NF}; /perf:mesh '"$mesh"' /{print c,$NF}' $input_clean > mesh_$mesh.data
    fi
done

for balance in $BALANCE_EFF; do
    if [[ ! -e "$balance.data" ]]; then
           echo "Generating $balance.data"
           awk '/Simulation cycle/{c=$NF}; /perf:balance '"$balance"' /{print c,$(NF-1)}' $input_clean > balance_$balance.data
       fi
done
for balance in $BALANCE_MAX; do
    if [[ ! -e "$balance.data" ]]; then
           echo "Generating $balance.data"
           awk '/Simulation cycle/{c=$NF}; /perf:balance '"$balance"' /{print c,$NF}' $input_clean > balance_$balance.data
       fi
done

# ==============================
# Generate plots from data files
# ==============================


$bindir/_plot-perf.py
