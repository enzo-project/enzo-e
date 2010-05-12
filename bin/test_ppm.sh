#!/bin/tcsh -f

set suffix = 

setenv LD_LIBRARY_PATH $HOME/lib:/opt/ibmcmp/xlf/13.1/lib64

foreach opt (yes)
   set n = 3
   while ($n <= 128)

     if ($n < 10) then
        set cycles = 10
     else
        set cycles = 10
     endif
     ./test_ppm implosion3 $n $cycles 0 >> out.test_ppm${suffix}

     @ n = $n + 1
   end
end

