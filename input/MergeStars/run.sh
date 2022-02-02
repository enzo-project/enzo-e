#!/bin/bash

rm -r run-*

if [ ! -f particles.dat ]
then
    echo "Generating ICs"
    
    python ics.py -r 4.0e16 -m 2.0e36 -c 1.5e17 1.5e17 1.5e17 -d 2.0e8 2.0e8 2.0e8 -i 1.0e8 -n 1000
fi

echo "Running the Merge Stars Test"

$CHARM_HOME/bin/charmrun +p1 ++local $ENZO_E_HOME/bin/enzo-e merge_stars_test-1.in > merge_stars_test-1.out 2>&1

$CHARM_HOME/bin/charmrun +p8 ++local $ENZO_E_HOME/bin/enzo-e merge_stars_test-8.in > merge_stars_test-8.out 2>&1

$CHARM_HOME/bin/charmrun +p16 ++local $ENZO_E_HOME/bin/enzo-e merge_stars_test-16.in > merge_stars_test-16.out 2>&1 


mpirun -np 16 python images.py -i Dir_Merge-Stars-1 -o image-1 > images-1.out 2>&1
mpirun -np 16 python images.py -i Dir_Merge-Stars-8 -o image-8 > images-8.out 2>&1
mpirun -np 16 python images.py -i Dir_Merge-Stars-16 -o image-16 > images-16.out 2>&1

mkdir run-1
mv image-1_* run-1/
mkdir run-8
mv image-8_* run-8/
mkdir run-16
mv image-16_* run-16/

mpirun -np 16 python mass_momentum_conservation.py -i Dir_Merge-Stars-1 -o mmc-1.png > mmc-1.out 2>&1

mpirun -np 16 python mass_momentum_conservation.py -i Dir_Merge-Stars-8 -o mmc-8.png > mmc-8.out 2>&1

mpirun -np 16 python mass_momentum_conservation.py -i Dir_Merge-Stars-16 -o mmc-16.png > mmc-16.out 2>&1

mv mmc-1.png run-1/
mv mmc-8.png run-8/
mv mmc-16.png run-16/

rm -r Dir*
