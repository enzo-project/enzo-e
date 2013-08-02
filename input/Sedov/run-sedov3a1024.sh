#!/bin/bash
#PBS -q normal
#PBS -l nodes=64:ppn=16:native
#PBS -l walltime=0:30:00
#PBS -N sedov3a1024
#PBS -o out.stdout
#PBS -e out.stderr
#PBS -M jobordner@ucsd.edu
#PBS -m abe
#PBS -V
# Start of user commands - comments start with a hash sign (#)

P=1024
T=3a

source $HOME/Cello/cello-src/input/Sedov/include.sh
