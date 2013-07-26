#!/bin/bash
#PBS -q normal
#PBS -l nodes=1:ppn=16:native
#PBS -l walltime=1:00:00
#PBS -N sedov2a0001
#PBS -o out.stdout
#PBS -e out.stderr
#PBS -M jobordner@ucsd.edu
#PBS -m abe
#PBS -V
# Start of user commands - comments start with a hash sign (#)

cd $HOME/Cello/cello-src/input/Sedov

P=0001
charmrun=$HOME/Charm/charm/bin/charmrun

$charmrun ++mpiexec +p$P bin/charm/enzo-p sedov2a$P.in >& out.sedov2a$P


