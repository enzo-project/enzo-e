#!/bin/bash
#PBS -q normal
#PBS -l nodes=4:ppn=16:native
#PBS -l walltime=1:00:00
#PBS -N sedov3a0064
#PBS -o out.stdout
#PBS -e out.stderr
#PBS -M jobordner@ucsd.edu
#PBS -m abe
#PBS -V
# Start of user commands - comments start with a hash sign (#)

cd $HOME/Cello/cello-src/

P=0064

charmrun ++mpiexec +p$P bin/charm/enzo-p input/sedov3a$P.in >& out.sedov3a$P


