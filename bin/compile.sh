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

printf "cleaning all..."
foreach type ($types)
   scons arch=$arch type=$type -c >& /dev/null
end
printf "\n\n"

foreach type ($types)

   set platform = $arch-$type

   set time = `date +"%Y-%m-%d %H:%M:%S"`
   printf "$time "
   printf "%14s" "${platform}: "

   # COMPILE

   scons arch=$arch type=$type -c >& /dev/null

   set start_hour = `date +"%H"`
   set start_min  = `date +"%M"`
   set start_sec  = `date +"%S"`

   scons -k -j$p arch=$arch type=$type >& out.scons.$platform

   @ hours = `date +"%H"` - $start_hour
   @ mins  = `date +"%M"` - $start_min
   @ secs  = `date +"%S"` - $start_sec

   @ total_secs = $secs + 60 * ( $mins + 60 * $hours )

   printf "%4d s  " $total_secs

   # count crashes

   grep FAIL test/$type-*unit |sort > fail.$platform
   grep pass test/$type-*unit |sort > pass.$platform
   set p = `cat pass.$platform | wc -l`
   set f = `cat fail.$platform | wc -l`

   printf "FAIL: %d  Pass: %d  " $f $p

   # check if any tests didn't finish

   set crash = `grep "UNIT TEST" test/$type-*unit | sed 's/BEGIN/END/' | uniq -u | wc -l`

   printf "crash: %d\n" $crash
   if ($crash != 0) then
      grep "UNIT TEST" test/$type-*unit | sed 's/BEGIN/END/' | uniq -u | sed 's/:/ /' | awk '{print "   ", $1}'
   endif


end

set time = `date +"%Y-%m-%d %H:%M:%S"`

echo "$time"

grep Segmentation out.*

echo


