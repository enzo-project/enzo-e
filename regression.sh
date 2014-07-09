#!/bin/tcsh -f

set target = test

set H0 = `date +"%H"`
set M0 = `date +"%M"`
set S0 = `date +"%S"`

# ./compile.sh clean

set H = `date +"%H"`
set M = `date +"%M"`
set S = `date +"%S"`

@ t = ($S - $S0) + 60 * ( ( $M - $M0) + 60 * ( $H - $H0) )
echo "END  time = $t s"


   set arch = $CELLO_ARCH
   set prec = $CELLO_PREC
   
   echo "BEGIN arch = $arch  prec = $prec"

   set H0 = `date +"%H"`
   set M0 = `date +"%M"`
   set S0 = `date +"%S"`

   ./compile.sh $target

   cp test/out.scons out.scons.$arch-$prec

   set H = `date +"%H"`
   set M = `date +"%M"`
   set S = `date +"%S"`

   @ t = ($S - $S0) + 60 * ( ( $M - $M0) + 60 * ( $H - $H0) )

   echo "END  arch = $arch  prec = $prec  time = $t s"
# end
# end
