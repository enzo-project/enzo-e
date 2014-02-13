#!/bin/bash

function usec {

   awk '/Simulation cycle/{c = $5;}; /Performance '"$2"' time-usec/{print c,$6}' < $1 > perf/usec-$2.data
}
# 0 00005.61 Performance initial time-usec 7202150
# 0 00005.61 Performance adapt time-usec 0
# 0 00005.61 Performance refresh time-usec 977939
# 0 00005.61 Performance compute time-usec 0
# 0 00005.61 Performance output time-usec 16077497
# 0 00005.61 Performance stopping time-usec 347888
mkdir perf

usec $1 simulation
usec $1 initial
usec $1 cycle
usec $1 adapt
usec $1 refresh
usec $1 compute
usec $1 output
usec $1 stopping

cd perf
ch-plot-usec.py

