#!/bin/tcsh -f

set ARCH = $CELLO_ARCH
set TYPE = $CELLO_TYPE
set PREC = $CELLO_PREC

set proc = 8

# set TYPE = (mpi charm serial)
# set PREC = (single double)

if ($#argv >= 1) then
   if ($argv[1] == "clean") then
      set d = `date +"%Y-%m-%d %H:%M:%S"`
      printf "$d %-14s" "cleaning..."
      foreach type (serial mpi charm)
         scons arch=$ARCH type=$type -c >& /dev/null
         rm -rf test/$type >& /dev/null
         rm -rf bin/$type >& /dev/null
         rm -rf lib/$type >& /dev/null
      end
      rm -rf include >& /dev/null
      rm -f test/COMPILING
      rm -f input/*.in.out >& /dev/null
      rm -rf build
      printf "done\n"
      exit
   endif
endif

echo "ARCH = $ARCH"
echo "TYPE = $TYPE"
echo "PREC = $PREC"

set d = `date +"%Y-%m-%d %H:%M:%S"`
echo "$d BEGIN"

foreach arch ($ARCH)
foreach type ($TYPE)
foreach prec ($PREC)

   rm -f "test/*/running.$arch.$prec"

   set configure = $arch-$type-$prec


   printf "$type" > test/COMPILING

   # clean
  scons arch=$arch type=$type prec=$prec -c >& /dev/null

   # COMPILE

   set d = `date +"%Y-%m-%d %H:%M:%S"`

   printf "$d %-14s %-14s" "$arch $type $prec" "compiling..."

   set dir = test/$type

   if (! -d $dir) mkdir $dir

   touch "$dir/running.$arch.$prec"

   # compile code with -j $proc processors
   scons arch=$arch type=$type prec=$prec -j$proc -k enzo-p >& $dir/out.scons

   # need to finish running tests with -j 1
   scons arch=$arch type=$type prec=$prec -k  >>& $dir/out.scons
   rm -f "$dir/running.$arch.$prec"
  
   printf "done\n"

   # count crashes

   cat $dir/*unit |grep FAIL      | grep "0/" | sort > $dir/fail.$configure
   cat $dir/*unit |grep incomplete| grep "0/" | sort > $dir/incomplete.$configure
   cat $dir/*unit |grep pass      | grep "0/" | sort > $dir/pass.$configure
   set f = `cat $dir/fail.$configure | wc -l`
   set i = `cat $dir/incomplete.$configure | wc -l`
   set p = `cat $dir/pass.$configure | wc -l`

   set d = `date +"%Y-%m-%d %H:%M:%S"`

   set line = "$d ${configure} FAIL: $f Incomplete: $i Pass: $p "

   set crash = `grep "UNIT TEST" $dir/*unit | sed 's/BEGIN/END/' | uniq -u | wc -l`

   if ($crash != 0) then
      line = "$line CRASH: $crash"
      grep "UNIT TEST" $dir/*unit \
       | sed 's/BEGIN/END/ ; s/:/ /' \
       | uniq -u \
       | awk '{print "   ", $1}'
   endif

   echo $line
   echo $line >> compile.log

   rm -f test/COMPILING

end
end
end

set d = `date +"%Y-%m-%d %H:%M:%S"`
echo "$d END"



