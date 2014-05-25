#!/bin/tcsh -f

foreach file (*.png)
   echo $file
#   convert -resize 1600 $file temp.png
#   convert -crop 1600x900+0+350 temp.png $file:r-crop.png
end

png2swf *mesh*.png -r 20 -o mesh.swf
png2swf *de*.png -r 20 -o de.swf
png2swf *te*.png -r 20 -o te.swf


