#!/bin/tcsh

if ($#argv != 3) then
   echo "\n  Usage: $0 <arch> <type> <prec>\n"
   exit
else 
   set arch = $argv[1]
   set type = $argv[2]
   set prec = $argv[3]
endif

echo
echo "arch = $arch"
echo "type = $type"
echo "prec = $prec"
echo

if ($arch == "triton-gnu") then
    module remove intel pgi
    module add gnu
    module add hdf5 mpich2_mx
else if ($arch == "triton-pgi") then
    module remove gnu intel
    module add pgi
    module add hdf5 mpich2_mx
else if ($arch == "triton-intel") then
    module remove pgi gnu
    module add intel
    module add hdf5 mpich2_mx
endif

set procs = 1

rm -f "test/*/running.$arch.$prec"

set d = `date +"%Y-%m-%d %H:%M:%S"`
echo "$d BEGIN"

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
     rm -f input/*.in.out >& /dev/null
     rm -rf build
     printf "done\n"

   else

      printf "$type" > test/COMPILING

      # COMPILE

      set d = `date +"%Y-%m-%d %H:%M:%S"`

      printf "$d %-14s %-14s" "${platform}" "compiling..."

      if (! -d $dir) mkdir $dir

      touch "$dir/running.$arch.$prec"

      scons arch=$arch type=$type -k -j8 bin/$type/enzo-p >>& $dir/out.scons
      scons arch=$arch type=$type -k -j$procs             >>& $dir/out.scons
      scons arch=$arch type=$type -k -j$procs             >>& $dir/out.scons
      rm -f "$dir/running.$arch.$prec"
  
      printf "done\n"

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
            grep "UNIT TEST" $dir/*unit \
             | sed 's/BEGIN/END/ ; s/:/ /' \
             | uniq -u \
             | awk '{print "   ", $1}'
         endif
      else 
         printf "\n"
      endif

      rm -f test/COMPILING

   endif


set d = `date +"%Y-%m-%d %H:%M:%S"`
echo "$d END"



