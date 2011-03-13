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
  set types = (serial mpi charm )
endif

set prec = $CELLO_PREC

echo
echo "arch  = $arch"
echo "types = ( $types )"
echo "prec  = $prec"
echo
set procs = 1

rm running.*.*.*

foreach type ($types)

   set platform = $arch-$type

   set d = `date +"%Y-%m-%d %H:%M:%S"`

   printf "$d %-14s %-14s" "${platform}" "cleaning..."
   scons arch=$arch type=$type -c >& /dev/null
   rm -f test/$type-*unit >& /dev/null
   rm -f bin-$type/* >& /dev/null
   printf "done\n"


   # COMPILE

   set d = `date +"%Y-%m-%d %H:%M:%S"`

   printf "$d %-14s %-14s" "${platform}" "compiling..."

   touch "running.$arch.$type.$prec"

   set t = `(time scons arch=$arch type=$type -k -j$procs >& out.scons.$platform)`
   rm -f "running.$arch.$type.$prec"
  
   set secs = `echo $t | awk '{print $3}'`

   printf "done (%4s s)\n" $secs

   # count crashes

   cat test/$type-*unit | grep FAIL | grep "0/" | sort > fail.$platform
   cat test/$type-*unit | grep pass | grep "0/" | sort > pass.$platform
   cat test/$type-*unit | grep incomplete | grep "0/" | sort > incomplete.$platform
   set p = `cat pass.$platform | wc -l`
   set f = `cat fail.$platform | wc -l`
   set i = `cat incomplete.$platform | wc -l`

   set d = `date +"%Y-%m-%d %H:%M:%S"`

   printf "$d %-14s " ${platform} 

   printf "FAIL: $f "

   printf "Pass: $p "

   printf "Incomplete: $i "

   # check if any tests didn't finish

   set crash = `grep "UNIT TEST" test/$type-*unit | sed 's/BEGIN/END/' | uniq -u | wc -l`

   if ($crash != 0) then
      printf "CRASH: $crash"
      if ($crash != 0) then
         grep "UNIT TEST" test/$type-*unit | sed 's/BEGIN/END/' | uniq -u | sed 's/:/ /' | awk '{print "   ", $1}'
      endif
   endif

   printf "\n"


end

set d = `date +"%Y-%m-%d %H:%M:%S"`
echo "$d"

grep Segmentation out.*

echo


