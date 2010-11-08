#!/bin/csh -f

rm -f out.scons.* fail.scons.*

if ($#argv >= 2) then
  set arch = ($argv[1])
  set types = ($argv[2-])
else if ($#argv == 1) then
  set arch = $CELLO_ARCH
  set types = ($argv)
else
  set arch = $CELLO_ARCH
  set types = (charm serial ampi mpi)
endif

echo
echo "arch  = $arch"
echo "types = ( $types )"
echo
set p = 1

# clear

printf "cleaning..."
foreach type ($types)
   setenv CELLO_TYPE $type
   scons -c >& /dev/null
end
printf "\n\n"

foreach type ($types)

   set platform = $arch-$type

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

set time = `date +"%Y-%m-%d %H:%M:%S"`

echo "$time"

grep Segmentation out.*

echo


