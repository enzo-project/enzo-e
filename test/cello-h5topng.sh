#!/bin/bash

test=$1
cycle=$2

min=0.0
max=1.0

echo "$test-*-$cycle.h5"
file_format="$test-*-$cycle.h5"
for file in $file_format;
do
    echo $file
    for block in `h5ls $file | awk '{print $1}'`;
    do
	h5topng -m $min -M $max -o $test-$cycle-$block.png ${file}:"$block/field density"
    done
done
