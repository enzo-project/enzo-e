#!/bin/tcsh -f

set ARCH = $CELLO_ARCH
set PREC = $CELLO_PREC

set log = "log.compile"
set proc = 8

set target = "install-bin"

if ($#argv >= 1) then
   if ($argv[1] == "clean") then
      set d = `date +"%Y-%m-%d %H:%M:%S"`
      printf "$d %-14s" "cleaning..."
      foreach prec (single double)
         python scons.py arch=$ARCH -c >& /dev/null
         rm -rf bin >& /dev/null
         rm -rf lib >& /dev/null
      end
      rm -rf include >& /dev/null
      rm -f test/COMPILING
      rm -f input/*.in.out >& /dev/null
      rm -rf build
      printf "done\n"
      exit
   else if ($argv[1] == "compile") then
     set target = install-bin
   else if ($argv[1] == "test") then
     set target = ""
   else
     # assume enzo-p
     set target = $argv[1]
   endif
endif

echo "ARCH = $ARCH"
echo "PREC = $PREC"
echo "target = $target"

set d = `date +"%Y-%m-%d %H:%M:%S"`
echo "$d BEGIN"

foreach arch ($ARCH)
foreach prec ($PREC)

    echo "arch = $arch  prec = $prec"

   rm -f "test/*/running.$arch.$prec"

   set configure = $arch-$prec
   set configure_print = `printf "%s %s %s" $arch $prec`


   printf "$configure" > test/COMPILING

   # make output directory for compilation and tests

   set dir = test

   # COMPILE

   set d = `date +"%Y-%m-%d %H:%M:%S"`

   printf "$d %-14s %-14s" "$arch $prec" "compiling..."
   printf "$d %-14s %-14s" "$arch $prec" "compiling..." >> $log

   touch "$dir/running.$arch.$prec"

   setenv CELLO_ARCH $arch
   setenv CELLO_PREC $prec

   rm -f bin/enzo-p

   python scons.py install-inc    >&  $dir/out.scons
   python scons.py -j $proc -Q $target |& tee $dir/out.scons

   util/parse_error.awk   < $dir/out.scons >  errors.org
   util/parse_warning.awk < $dir/out.scons >> errors.org

   if (-e bin/enzo-p) then
      echo "Success"
   else
      echo "FAIL"
   endif

   rm -f "$dir/running.$arch.$prec"

   printf "done\n"
   printf "done\n" >> $log

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

      printf "%s %s %-12s %-6s %-6s %s %-2s %s %-2s %s %-4s %s %-2s\n" $line
      printf "%s %s %-12s %-6s %-6s %s %-2s %s %-2s %s %-4s %s %-2s\n" $line >> $log

      foreach test ($dir/*unit)
        set test_begin = `grep "UNIT TEST BEGIN" $test | wc -l`
        set test_end   = `grep "UNIT TEST END" $test | wc -l`

        @ crash = $test_begin - $test_end

        if ($crash != 0) then
           set line = "   CRASH: $test\n"
           printf "$line"
           printf "$line" >> $log
        endif
      end


   endif

   if ($CELLO_ARCH == "ncsa-bw") then
       echo "Relinking with static libpng15.a..."
      /u/sciteam/bordner/Charm/charm/bin/charmc -language charm++ -tracemode projections -o build/charm/Enzo/enzo-p -g -g build/charm/Enzo/enzo-p.o build/charm/Cello/main_enzo.o -Llib/charm -L/opt/cray/hdf5/default/cray/74/lib -lcharm -lenzo -lsimulation -lproblem -lcomm -lmesh -lfield -lio -ldisk -lmemory -lparameters -lerror -lmonitor -lparallel -lperformance -ltest -lcello -lexternal -lhdf5 -lz /u/sciteam/bordner/lib/libpng15.a -lgfortran

       mv build/charm/Enzo/enzo-p bin/charm

   endif
   rm -f test/COMPILING

end
end

set d = `date +"%Y-%m-%d %H:%M:%S"`
echo "$d END"



