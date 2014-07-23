#!/bin/bash

#foreach file (*.png)
#   echo $file
#   convert -resize 1600 $file temp.png
#   convert -crop 1600x900+0+350 temp.png $file:r-crop.png
#end

out=$1
shift
png2swf $@ -r 20 -o $out.swf


