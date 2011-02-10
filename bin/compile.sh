#!/bin/tcsh -f

rm -f out.scons.* fail.scons.*

if ($#argv >= 2) then
  set arch = ($argv[1])
  set types = ($argv[2-])
else if ($#argv == 1) then
  set arch = $CELLO_ARCH
  set types = ($argv)
else
  set arch = $CELLO_ARCH
  set types = (serial mpi ampi charm )
endif

echo
echo "arch  = $arch"
echo "types = ( $types )"
echo
set procs = 4

# clear

set d = `date +"%Y-%m-%d %H:%M:%S"`
printf "$d cleaning all..."
foreach type ($types)
   scons arch=$arch type=$type -c >& /dev/null
end
printf "\n\n"

foreach type ($types)

   set platform = $arch-$type

   set d = `date +"%Y-%m-%d %H:%M:%S"`
   printf "$d %14s" "${platform}: "

   # COMPILE

   scons arch=$arch type=$type -c >& /dev/null

   set t = `(time scons arch=$arch type=$type -k -j$procs >& out.scons.$platform)`
  
   set secs = `echo $t | awk '{print $3}'`

   printf "%4s s  " $secs

   # count crashes

   cat test/$type-*unit | grep FAIL | grep "0/" | sort > fail.$platform
   cat test/$type-*unit | grep pass | grep "0/" | sort > pass.$platform
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

set d = `date +"%Y-%m-%d %H:%M:%S"`
echo "$d"

grep Segmentation out.*

echo


