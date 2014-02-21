#!/bin/bash

if ("x$2" == "x"); then
    echo "Usage: $0 <enzo-p output>"
    exit 1
fi
function usec {
   awk '/Simulation cycle/{c = $5;}; /Performance '"$2"' time-usec/{print c,$6}' < $1 > perf/time-$2.data
}

function mem {
   awk '/Simulation cycle/{c = $5;}; /Performance cycle '"$2" '/{print c,$6}' < $1 > perf/mem-$2.data
}

mkdir perf

usec $1 simulation
usec $1 initial
usec $1 cycle
usec $1 adapt
usec $1 refresh
usec $1 compute
usec $1 output
usec $1 stopping

# mem $1 "bytes-curr"
cd perf
$HOME/Cello/cello-src/util/ch-perf.py

