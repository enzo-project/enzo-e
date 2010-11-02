#!/bin/csh -f

rm -f out.scons.* fail.scons.*

setenv CELLO_ARCH linux

if ($#argv >= 1) then
  set types = ($argv)
else
  set types = (charm serial ampi mpi)
endif

printf "\n[ $types ]\n\n"
set p = 1

# clear

printf "cleaning..."
foreach type ($types)
   setenv CELLO_TYPE $type
   scons -c >& /dev/null
end
printf "\n\n"

foreach type ($types)

   setenv CELLO_TYPE $type

   set platform = ${CELLO_ARCH}-${CELLO_TYPE}

   set time = `date +"%Y-%m-%d %H:%M:%S"`
   printf "$time "
   printf "%10s" "${platform}: "
   # compile
   scons -k -j$p >& out.scons.$platform
   # count fails
   grep FAIL test/$type-*unit |sort > fail.$platform
   grep pass test/$type-*unit |sort > pass.$platform
   set p = `cat pass.$platform | wc -l`
   set f = `cat fail.$platform | wc -l`
   printf "FAIL: %d  Pass: %d\n" $f $p
   # check if any tests didn't finish
   grep "UNIT TEST" test/$type-*unit | sed 's/BEGIN/END/' | uniq -u

end

echo

grep Segmentation out.*

echo


