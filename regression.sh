#!/bin/csh -f

set target = test

./compile.sh clean

foreach type (charm mpi)
foreach arch (gnu)
foreach prec (double)

   
   echo "BEGIN type = $type  arch = $arch  prec = $prec"
   setenv CELLO_ARCH linux-$arch
   setenv CELLO_TYPE $type
   setenv CELLO_PREC $prec

   set H0 = `date +"%H"`
   set M0 = `date +"%M"`
   set S0 = `date +"%S"`

   ./compile.sh $target

   cp test/$type/out.scons out.scons.$type-$arch-$prec

   set H = `date +"%H"`
   set M = `date +"%M"`
   set S = `date +"%S"`

   echo "$H0 $M0 $S0"
   echo "$H $M $S"

   @ t = ($S - $S0) + 60 * ( ( $M - $M0) + 60 * ( $H - $H0) )

   echo "END  type = $type  arch = $arch  prec = $prec  time = $t"
end
end
end
