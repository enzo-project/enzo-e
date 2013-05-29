#!/bin/tcsh -f

set target = test

set H0 = `date +"%H"`
set M0 = `date +"%M"`
set S0 = `date +"%S"`

./compile.sh clean

set H = `date +"%H"`
set M = `date +"%M"`
set S = `date +"%S"`

@ t = ($S - $S0) + 60 * ( ( $M - $M0) + 60 * ( $H - $H0) )
echo "END  time = $t s"


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

   @ t = ($S - $S0) + 60 * ( ( $M - $M0) + 60 * ( $H - $H0) )

   echo "END  type = $type  arch = $arch  prec = $prec  time = $t s"
end
# end
# end
