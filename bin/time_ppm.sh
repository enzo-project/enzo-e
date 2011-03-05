#!/bin/tcsh -f

set suffix = $1

set nmin   = 2
set nmax   = 100
set cycles = 10

# setenv LD_LIBRARY_PATH $HOME/lib:/opt/ibmcmp/xlf/13.1/lib64

foreach opt (yes)
   set n = $nmin
   while ($n <= $nmax)

      echo $n
     ./test_ppm${suffix} ppm-implosion3 $n $cycles 0 > test_ppm.out
      awk '/Time proc /{print $4}' test_ppm.out >> out.time-ppm${suffix}
      awk '/GFlop rate/{print $4}' test_ppm.out >> out.gflop-ppm${suffix}
      awk '/Flop count/{print $4}' test_ppm.out >> out.flops-ppm${suffix}

     @ n = $n + 1
   end
end

