#!/bin/tcsh -f

set n = 004
while ($n <= 200)
   m4 -DN=$n input/test_EnzoMethodPpm.in.m4 > input/test_EnzoMethodPpm.in

   set mflop_rate = `bin/test_EnzoMethodPpm | awk '/MFlop/{print $4}'`
   printf "$n $mflop_rate\n"

   @ n = $n + 1
end

     
