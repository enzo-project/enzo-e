#!/bin/tcsh -f

set suffix = $1

set nmin   = 4
set nmax   = 64
set cycles = 10

# setenv LD_LIBRARY_PATH $HOME/lib:/opt/ibmcmp/xlf/13.1/lib64

foreach opt (yes)
   set n = $nmin
   while ($n <= $nmax)

     ./test_ppml ppml-implosion3 $n $cycles 0 | grep time >> out.test_ppml${suffix}

     @ n = $n + 1
   end
end

