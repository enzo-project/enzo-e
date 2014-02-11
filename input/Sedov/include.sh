#----------------------------------------------------------------------
# INPUT FROM CALLING SCRIPT
#
#    P=%04d     number of processors
#    T=[2a|3a]  type: 2D amr or unigrid
#    H=["sdsc-gordon"|"sdsc-bw"]
#
#----------------------------------------------------------------------

cello=$HOME/Cello/cello-src
charm=$HOME/Charm/charm

cd $cello/input/Sedov

enzorun=$cello/bin/enzo-p
charmrun=$charm/bin/charmrun
input=sedov$T$P.in
output=out.sedov$T$P

if [ $H == "sdsc-gordon" ]; then

   $charmrun ++mpiexec +p$P $enzorun $input >& $output

elif [ $H == "ncsa-bw" ]; then

    aprun -n 64 $enzorun $input >& $output

fi



