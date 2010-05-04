#!/bin/tcsh -f

set suffix = 


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

