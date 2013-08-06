#----------------------------------------------------------------------
# INPUT FROM CALLING SCRIPT
#
#    P=%04d     number of processors
#    T=[2a|3a]  type: 2D amr or unigrid
#
#----------------------------------------------------------------------
# FILES
#
#  CELLO: e.g. cello-src.infiniband
#  CHARM: e.g. 650/gnu/infiniband/charm-6.5.0
#
#----------------------------------------------------------------------

cello=$HOME/Cello/cello-src
charm=$HOME/Charm/charm

cd $cello/input/Sedov

enzorun=$cello/bin/charm/enzo-p
charmrun=$charm/bin/charmrun
input=sedov$T$P.in
output=out.sedov$T$P

$charmrun ++mpiexec +p$P $enzorun $input >& $output
# $charmrun +p$P $enzorun $input >& $output


