#!/bin/bash
#PBS -l nodes=1:ppn=32:xe
#PBS -l walltime=00:10:00
#PBS -N testjob

P=16
H="ncsa-bw"

source $HOME/Cello/cello-src/input/Sedov/include.sh
