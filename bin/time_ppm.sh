#!/bin/tcsh -f

set suffix = $1

set nmin   = 4
set nmax   = 100
set cycles = 10

# setenv LD_LIBRARY_PATH $HOME/lib:/opt/ibmcmp/xlf/13.1/lib64

foreach opt (yes)
   set n = $nmin
   while ($n <= $nmax)

     ./test_ppm ppm-implosion3 $n $cycles 0 | grep time  >> out.time-ppm.${suffix}
     ./test_ppm ppm-implosion3 $n $cycles 0 | grep GFlop >> out.gflop-ppm.${suffix}

     @ n = $n + 1
   end
end

