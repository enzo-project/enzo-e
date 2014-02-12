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

P0=`printf "%04d" $P`
enzorun=$cello/bin/enzo-p
charmrun=$charm/bin/charmrun
input=sedov$T$P0.in
output=out.sedov$T$P0



if [ $H == "sdsc-gordon" ]; then

   $charmrun ++mpiexec +p$P $enzorun $input >& $output

elif [ $H == "ncsa-bw" ]; then


    . /opt/modules/default/init/bash
    module swap PrgEnv-cray PrgEnv-gnu

   aprun -n $P -d 2 $enzorun $input >& $output

fi



