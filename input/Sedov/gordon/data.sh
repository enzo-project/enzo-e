#!/bin/bash

rm -f *.data

for file in out.sedov*
do
   P=`awk '/Simulation processors/{print $6}' $file`
   echo $P >> procs.data
   awk '/Simulation cycle 0000/{print $2}' $file >> time-startup.data
   awk '/Simulation cycle 0010/{print $2}' $file >> time-total.data
   bytes_n=`cat $file | awk '/Performance cycle bytes-highest/{print $6}' $file | tail -1`
   echo "scale=10; $bytes_n / $P" | bc >> bytes-highest.data
done

python mem-total.py
python time-total.py


