#!/bin/csh -f

set target = test

./compile.sh clean

foreach type (charm mpi)

   set arch = $CELLO_ARCH
   set prec = $CELLO_PREC

   
   echo "BEGIN type = $type  arch = $arch  prec = $prec"
   setenv CELLO_TYPE $type

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

   echo "END  type = $type  arch = $arch  prec = $prec  time = $t s"
end
# end
# end
