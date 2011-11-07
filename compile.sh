#!/bin/tcsh -f

set ARCH = $CELLO_ARCH
set TYPE = $CELLO_TYPE
set PREC = $CELLO_PREC

set proc = 8

# set TYPE = (mpi charm serial)
# set PREC = (single double)

# set default target
set target = "compile"

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
   else if ($argv[1] == "enzo-p") then
     set target = bin/$type/enzo-p
   else if ($argv[1] == "compile") then
     set target = install-bin
   else if ($argv[1] == "test") then
     set target = ""
   else
     echo "Usage: $0 [clean|compile|test|enzo-p]"
     exit(1)
   endif
endif

echo "ARCH = $ARCH"
echo "TYPE = $TYPE"
echo "PREC = $PREC"
echo "target = $target"

set d = `date +"%Y-%m-%d %H:%M:%S"`
echo "$d BEGIN"

foreach arch ($ARCH)
foreach type ($TYPE)
foreach prec ($PREC)

   rm -f "test/*/running.$arch.$prec"

   set configure = $arch-$type-$prec
   set configure_print = `printf "%s %s %s" $arch $type $prec`


   printf "$type" > test/COMPILING

   # clean
#  scons arch=$arch type=$type prec=$prec -c >& /dev/null

   # make output directory for compilation and tests

   set dir = test/$type
   if (! -d $dir) mkdir $dir

   # COMPILE

   set d = `date +"%Y-%m-%d %H:%M:%S"`

   printf "$d %-14s %-14s" "$arch $type $prec" "compiling..."

   touch "$dir/running.$arch.$prec"

   scons arch=$arch type=$type prec=$prec -k -j$proc install-bin >& $dir/out.scons
   scons arch=$arch type=$type prec=$prec -k         $target    >>& $dir/out.scons

   rm -f "$dir/running.$arch.$prec"

   printf "done\n"

   # TESTS

   if ($target == "") then
  
      # count crashes

      cat $dir/*unit |grep FAIL      | grep "0/" | sort > $dir/fail.$configure
      cat $dir/*unit |grep incomplete| grep "0/" | sort > $dir/incomplete.$configure
      cat $dir/*unit |grep pass      | grep "0/" | sort > $dir/pass.$configure
      set f = `cat $dir/fail.$configure | wc -l`
      set i = `cat $dir/incomplete.$configure | wc -l`
      set p = `cat $dir/pass.$configure | wc -l`

      set d = `date +"%Y-%m-%d %H:%M:%S"`

      set line = "$d ${configure_print} FAIL: $f Incomplete: $i Pass: $p "

      set crash = `grep "UNIT TEST" $dir/*unit | sed 's/BEGIN/END/' | uniq -u | wc -l`
      if ($crash != 0) then
         set line = "$line CRASH: $crash"
         grep "UNIT TEST" $dir/*unit \
          | sed 's/BEGIN/END/ ; s/:/ /' \
          | uniq -u \
          | awk '{print "   ", $1}'
      endif

      printf "%s %s %-12s %-6s %-6s %s %-2s %s %-2s %s %-4s %s %-2s\n" $line
      printf "%s %s %-12s %-6s %-6s %s %-2s %s %-2s %s %-4s %s %-2s\n" $line >> compile.log

   endif
   rm -f test/COMPILING

end
end
end

set d = `date +"%Y-%m-%d %H:%M:%S"`
echo "$d END"



