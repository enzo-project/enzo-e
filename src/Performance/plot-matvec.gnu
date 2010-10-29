set terminal png
set output "matvec.png"
set grid
set key bottom left

plot [][0:] './matvec.data' using ($1):($2) title "const" w l, \
            './matvec.data' using ($1):($3) title "symm" w l, \
            './matvec.data' using ($1):($4) title "full" w l, \
            './matvec.data' using ($1):($5) title "block" w l
