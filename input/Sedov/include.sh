#INPUT
# P=%04d     number of processors
# T=[2a|3a]  type: 2D amr or unigrid

cd $HOME/Cello/cello-src/input/Sedov

enzorun=$HOME/Cello/cello-src/bin/enzo-p
charmrun=$HOME/Charm/charm/bin/charmrun
input=input/sedov2a$P.in
output=out.sedov2a$P

$charmrun ++mpiexec +p$P bin/charm/enzo-p sedov2a$P.in >& out.sedov2a$P


