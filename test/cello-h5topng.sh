#!/bin/csh -f

set test = $1
set cycle = $2

set min = 0.0
set max = 1.0

echo "$test-*-$cycle.h5"
set file_format = "$test-*-$cycle.h5"
foreach file ($file_format)
   echo $file
   foreach block (`h5ls $file/patch_0 | awk '{print $1}'`)
      h5topng -m $min -M $max -o $test-$cycle-$block.png ${file}:patch_0/$block/density
   end
end
