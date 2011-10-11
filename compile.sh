#!/bin/tcsh -f

if ($#argv >= 1) then
  set types = ($argv)
else
   set types = (charm mpi serial)
endif

set arch = $CELLO_ARCH
set prec = $CELLO_PREC

echo
echo "arch  = $arch"
echo "types = ( $types )"
echo "prec  = $prec"
echo

set procs = 1

rm -f "test/*/running.$arch.$prec"

set d = `date +"%Y-%m-%d %H:%M:%S"`
echo "$d BEGIN"
foreach type ($types)

   set dir = test/$type
   set platform = $arch-$type

   set d = `date +"%Y-%m-%d %H:%M:%S"`

   if ($type == "clean") then


     printf "$d %-14s %-14s" "${platform}" "cleaning..."
     scons arch=$arch type=serial -c >& /dev/null
     scons arch=$arch type=mpi    -c >& /dev/null
     scons arch=$arch type=charm  -c >& /dev/null
     rm -f $dir/*unit >& /dev/null
     rm -f bin/$type/* >& /dev/null
     rm -f test/COMPILING
     rm -rf build
     printf "done\n"

   else

      printf "$type" > test/COMPILING

      # COMPILE

      set d = `date +"%Y-%m-%d %H:%M:%S"`

      printf "$d %-14s %-14s" "${platform}" "compiling..."

      if (! -d $dir) mkdir $dir

      touch "$dir/running.$arch.$prec"

      set t = `(time scons arch=$arch type=$type -k -j$procs >& $dir/out.scons)`
      rm -f "$dir/running.$arch.$prec"
  
      set secs = `echo $t | awk '{print $3}'`

      printf "done (%4s s)\n" $secs

      # count crashes

      cat $dir/*unit |grep FAIL      | grep "0/" | sort > $dir/fail.$platform
      cat $dir/*unit |grep incomplete| grep "0/" | sort > $dir/incomplete.$platform
      cat $dir/*unit |grep pass      | grep "0/" | sort > $dir/pass.$platform
      set f = `cat $dir/fail.$platform | wc -l`
      set i = `cat $dir/incomplete.$platform | wc -l`
      set p = `cat $dir/pass.$platform | wc -l`

      set d = `date +"%Y-%m-%d %H:%M:%S"`

      printf "$d %-14s " ${platform} 

      printf "FAIL: $f "

      printf "Incomplete: $i "

      printf "Pass: $p "

      # check if any tests didn't finish

      set crash = `grep "UNIT TEST" $dir/*unit | sed 's/BEGIN/END/' | uniq -u | wc -l`

      if ($crash != 0) then
         printf "CRASH: $crash"
         if ($crash != 0) then
            grep "UNIT TEST" $dir/*unit | sed 's/BEGIN/END/' | uniq -u | sed 's/:/ /' | awk '{print "   ", $1         }'
         endif
      endif

      printf "\n"

      rm -f test/COMPILING

   endif


end

set d = `date +"%Y-%m-%d %H:%M:%S"`
echo "$d END"



