set title  "Matvec Performance versus Matrix Storage Format"
set xlabel "N (problem size = N^3)"
set ylabel "MFlop/s"

set output  "plot-format.png"
set terminal png

plot 'matvec-krummhorn-O2-2.data' using ($1+2):($4) title 'Nonsymmetric' w l, \
     'matvec-krummhorn-O2-2.data' using ($1+2):($3) title 'Symmetric'  w l, \
     'matvec-krummhorn-O2-2.data' using ($1+2):($2) title 'Constant' w l

set output  "plot-format.eps"
set terminal postscript color eps font 18

plot 'matvec-krummhorn-O2-2.data' using ($1+2):($4) title 'Nonsymmetric' w l, \
     'matvec-krummhorn-O2-2.data' using ($1+2):($3) title 'Symmetric'  w l, \
     'matvec-krummhorn-O2-2.data' using ($1+2):($2) title 'Constant' w l

