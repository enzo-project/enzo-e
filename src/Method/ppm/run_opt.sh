#!/bin/tcsh

set n = 1
while ($n < 128)
  set c = `echo "2000000 / $n / $n" | bc`
  test_implosion.opt-no $n $c 0
  test_implosion.opt-yes $n $c 0
  @ n = $n + 1
end

