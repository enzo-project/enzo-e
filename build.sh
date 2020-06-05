#!/bin/bash
#
# Input environment variables
#
#   CELLO_ARCH
#   CELLO_PREC
#
# Output status files for ./build.sh test
#
#   test/STATUS
#   test/DATE
#   test/START
#   test/STOP

arch=$CELLO_ARCH
prec=$CELLO_PREC

python="python2"
# initialize time

H0=`date +"%H"`
M0=`date +"%M"`
S0=`date +"%S"`

log="log.build"

proc=4

# set default target

target="install-bin"
k_switch="-k"

if [ "$#" -ge 1 ]; then
   if [ "$1" == "clean" ]; then
       TMP=`mktemp -t build.sh.XXXXXX`
       trap "rm $TMP* 2>/dev/null" EXIT
       
       d=`date +"%Y-%m-%d %H:%M:%S"`
      printf "$d %-14s cleaning..."
      for prec in single double; do
         $python scons.py arch=$arch -c >& /dev/null
         rm -rf bin >& /dev/null
         rm -rf lib >& /dev/null
      done
      rm -rf include >& /dev/null
      rm -f input/*.in.out >& /dev/null
      rm -rf build build-*
      rm -rf test/*.h5
      rm -rf template_defs.def.h template_defs.decl.h
      rm -rf .sconf_temp/conftest_0.c .sconsign.dblite 
      rm -rf config.log warnings.org errors.org log.build out.scons.*
      rm -rf config/*.pyc
      rm -rf test/fail.* test/pass.* test/incomplete.*
      rm -rf test/*.test_log
      rm -rf scons-local-2.2.0/SCons/*.pyc scons-local-2.2.0/SCons/*/*.pyc
      rm -rf charmrun parameters.out checkpoint_ppm* output-stride*.h5
      rm -rf `find test -name "*.png"`
      rm -rf `find test -name "*.h5"`
      printf "done\n"
      rm -rf test/out.scons
   fi
      if [ "$1" == "reset" -o "$1" == "clean" ]; then
	cd test
	rm -f STATUS DATE START STOP ARCH PREC
       exit
   elif [ "$1" == "compile" ]; then
      target=install-bin
   elif [ "$1" == "test" ]; then
      ./build.sh
      target="test"
      proc=1
      k_switch="-k"
   elif [ "$1" == "help" ]; then
      echo
      echo "Usage: $0 [clean|compile|test]"
      echo
      echo "       $0 bin/enzo-p"
      echo
      echo "       $0 bin/test_Foo"
      echo
      echo "       $0 test/test_Foo.unit"
      exit
   else
      k_switch=""
      target="$1"
#      rm -f $target
	echo "Remove $target"
   fi
else
   # assume enzo-p
   k_switch=""
   target="bin/enzo-p"
fi

if [ $target == "bin/enzo-p" ]; then
   if [ -e $target ]; then
       echo "Saving existing bin/enzo-p to bin/enzo-p.prev"
       mv bin/enzo-p bin/enzo-p.prev
   fi
fi
    

echo "Compiling" > test/STATUS

date=`date +"%Y-%m-%d"`
start=`date +"%H:%M:%S"`
echo "$date $start BEGIN"

echo "BEGIN Enzo-P/Cello ${0}"
echo "arch=$arch"
echo "prec=$prec"
echo "target=$target"

rm -f "test/*/running.$arch.$prec"

configure=$arch-$prec
configure_print=`printf "%s %s %s" $arch $prec`

# make output directory for compilation and tests

dir=test

# COMPILE

d=`date +"%Y-%m-%d %H:%M:%S"`

printf "$date $start %-14s %-14s" "compiling..."
printf "$date $start %-14s %-14s" "compiling..." >> $log

touch "$dir/running.$arch.$prec"

CELLO_ARCH=$arch; export $CELLO_ARCH
CELLO_PREC=$prec; export $CELLO_PREC

if [ $target == "test" ]; then

    echo "$date"     > test/DATE
    echo "$start"    > test/START
    echo "$arch"     > test/ARCH
    echo "$prec"     > test/PREC
fi    


$python scons.py install-inc    &>  $dir/out.scons
$python scons.py $k_switch -j $proc -Q $target  2>&1 | tee $dir/out.scons

./tools/awk/error-org.awk   < $dir/out.scons >  errors.org
./tools/awk/warning-org.awk < $dir/out.scons >  warnings.org

if [ -e $target ]; then
   echo "Success"
else
   echo "FAIL"
fi

rm -f "$dir/running.$arch.$prec"

printf "done\n"
printf "done\n" >> $log

# TESTS

if [ $target == "test" ]; then

    rm -f              test/STOP

   # count failures, incompletes, and passes

   grep "^ FAIL"       $dir/*unit > $dir/fail.$configure
   grep "^ incomplete" $dir/*unit > $dir/incomplete.$configure
   grep "^ pass"       $dir/*unit > $dir/pass.$configure

   f=`wc -l < $dir/fail.$configure`
   i=`wc -l < $dir/incomplete.$configure`
   p=`wc -l < $dir/pass.$configure`

   stop=`date +"%H:%M:%S"`

   line="$stop ${configure_print} FAIL: $f Incomplete: $i Pass: $p "

   printf "%s %s %-12s %-6s %-6s %s %-2s %s %-2s %s %-4s %s %-2s\n" \
        $line | tee $log

   for test in $dir/*unit; do

      test_begin=`grep "UNIT TEST BEGIN" $test | wc -l`
      test_end=`grep "UNIT TEST END"   $test | wc -l`

      crash=$((test_begin - $test_end))

      if [ $crash != 0 ]; then
         line="   CRASH: $test\n"
         printf "$line"
         printf "$line" >> $log
      fi
   done

   echo "$stop" > test/STOP

fi

if [ x$CELLO_ARCH == "xncsa-bw" ]; then
    echo "Relinking with static libpng15.a..."
    build_dir="build"
   /u/sciteam/bordner/Charm/charm/bin/charmc -language charm++ -tracemode projections -o $build_dir/charm/Enzo/enzo-p -g -g $build_dir/charm/Enzo/enzo-p.o $build_dir/charm/Cello/main_enzo.o -Llib/charm -L/opt/cray/hdf5/default/cray/74/lib -lcharm -lenzo -lsimulation -lproblem -lcomm -lmesh -lfield -lio -ldisk -lmemory -lparameters -lerror -lmonitor -lparallel -lperformance -ltest -lcello -lexternal -lhdf5 -lz /u/sciteam/bordner/lib/libpng15.a -lgfortran

   mv $build_dir/charm/Enzo/enzo-p bin/charm

fi

cp test/out.scons out.scons.$arch-$prec

H1=`date +"%H"`
M1=`date +"%M"`
S1=`date +"%S"`

t=`echo "scale=2; (( $S1 - $S0 ) + 60 * ( ( $M1 - $M0 ) + 60 * ( $H1 - $H0) ))/60.0" | bc`

echo "END   Enzo-P/Cello ${0}: arch = $arch  prec = $prec  target = $target time = ${t} min"

d=`date "+%H:%M:%S"`

rm -f test/STATUS

if [ ! -e $target ]; then
   exit 1
fi

# check for failures 
if [ $target == "test" ]; then

    file_attempted=test/runs_attempted.$configure
    file_started=test/runs_started.$configure
    file_completed=test/runs_completed.$configure

    ls test/test_*.unit                   > $file_attempted
    grep -l "BEGIN" test/test_*.unit      > $file_started
    grep -l "END CELLO"  test/test_*.unit > $file_completed


    count_attempted=`cat $file_attempted | wc -l `
    count_started=`cat $file_started | wc -l`
    count_completed=`cat $file_completed | wc -l`
    
    echo "Test run summary"
    echo
    echo "   Test runs attempted: $count_attempted"
    echo "   Test runs started:   $count_started"
    if [ $count_attempted -gt $count_started ]; then
	echo "   --------------------"
	sort $file_attempted $file_started | uniq -u | awk '{print "   *** failed in startup: "$1}'
	echo "   --------------------"
    fi
    echo "   Test runs completed: $count_completed"
    if [ $count_started -gt $count_completed ]; then
	echo "   --------------------"
        sort $file_started $file_completed | uniq -u | awk '{print "   *** incomplete output: "$1}'
	echo "   --------------------"
    fi
    echo

    if [ $f -gt 0 ]; then
	echo "Exiting testing with failures:"
	echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
	cat $dir/fail.$configure
	echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
	exit_status=1
    else
	echo "Exiting testing with success"
	exit_status=0
    fi
fi

if [ $target = "test" ] && [ "$CELLO_PREC" = "double" ]; then
    # the vl+ct tests should be consolidated with the rest of the tests
    echo ""
    echo "--------------------"
    echo "Attempting to run VL+CT tests (only defined for double Precision)"
    ./test/run_vlct_test.sh
    result_code=$?
    if [ $result_code -gt 0 ]; then
        exit_status=1
    fi
fi;
echo "Done."
exit $exit_status
