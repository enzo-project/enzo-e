#!/bin/bash
#PBS -q normal
#PBS -l nodes=16:ppn=16:native
#PBS -l walltime=1:00:00
#PBS -N sedov3a0256
#PBS -o out.stdout
#PBS -e out.stderr
#PBS -M jobordner@ucsd.edu
#PBS -m abe
#PBS -V
# Start of user commands - comments start with a hash sign (#)

P=0256
T=3a

source $HOME/Cello/cello-src/input/Sedov/include.sh
