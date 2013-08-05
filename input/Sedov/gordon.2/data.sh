#!/bin/bash

rm -f *.data

for out in out.sedov3a????
do
   ip=`awk '/Simulation processors/{print $6}' $out`
   run=P$ip

   echo $ip >> procs.data

   awk '/Simulation cycle/ {print $2}'                    < $out > $run-wallclock.data
   awk '/simulation time-usec/{print $6/1000000/'"$ip"'}' < $out > $run-time-compute.data
   awk '/cycle time-usec/{print $6/1000000/'"$ip"'}'      < $out > $run-time-cycle.data
   awk '/initial time-usec/{print $6/1000000/'"$ip"'}'    < $out > $run-time-initial.data
   awk '/adapt time-usec/{print $6/1000000/'"$ip"'}'      < $out > $run-time-adapt.data
   awk '/refresh time-usec/{print $6/1000000/'"$ip"'}'    < $out > $run-time-refresh.data
   awk '/compute time-usec/{print $6/1000000/'"$ip"'}'    < $out > $run-time-compute.data
   awk '/output time-usec/{print $6/1000000/'"$ip"'}'     < $out > $run-time-output.data
   awk '/prepare time-usec/{print $6/1000000/'"$ip"'}'    < $out > $run-time-prepare.data
   awk '/simulation bytes-highest/{print $6}' < $out > $run-bytes-highest.data

done


