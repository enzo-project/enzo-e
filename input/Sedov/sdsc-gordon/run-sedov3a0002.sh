#!/bin/bash
#PBS -q normal
#PBS -l nodes=1:ppn=16:native
#PBS -l walltime=0:10:00
#PBS -N sedov3a0002
#PBS -o out.stdout
#PBS -e out.stderr
#PBS -M jobordner@ucsd.edu
#PBS -m abe
#PBS -V
# Start of user commands - comments start with a hash sign (#)

P=0002
T=3a
H="sdsc-gordon"

source $HOME/Cello/cello-src/input/Sedov/include.sh
