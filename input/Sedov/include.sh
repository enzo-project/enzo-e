#----------------------------------------------------------------------
# INPUT FROM CALLING SCRIPT
#
#    P=%04d     number of processors
#    H=["sdsc-gordon"|"sdsc-bw"]
#
#----------------------------------------------------------------------

P0=`printf "%04d" $P`

cello=$HOME/Cello/cello-src
charm=$HOME/Charm/charm
rundir=$cello/input/Sedov/run.$P0

rm -rf $rundir
mkdir $rundir
cd $rundir

enzorun=$cello/bin/enzo-p
charmrun=$charm/bin/charmrun
input=sedov.in
output=out.sedov

cp -r ../$input .
cp -r ../config .
cp -r ../config*-*.incl .
cp $enzorun .

if [ $H == "sdsc-gordon" ]; then

   $charmrun ++mpiexec +p$P ./enzo-p $input >& $output

elif [ $H == "ncsa-bw" ]; then

    . /opt/modules/default/init/bash
   source $HOME/bin/cmod.sh

   aprun -n $P -d 2 ./enzo-p $input >& $output

fi



