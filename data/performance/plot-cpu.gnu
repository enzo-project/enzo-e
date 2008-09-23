set title  "Matvec Performance versus Processor Type"
set xlabel "N (problem size = N^3)"
set ylabel "MFlop/s"

set output  "plot-cpu.png"
set terminal png

plot 'matvec-krummhorn-2.data' using ($1+2):($4) title 'Intel Core 2 Duo' w l, \
     'matvec-datastar-2.data'     using ($1+2):($3) title 'IBM Power4'  w l, \
     'matvec-padoan-2.data'       using ($1+2):($2) title 'AMD Opteron' w l

set output  "plot-cpu.eps"
set terminal postscript color eps font 18

plot 'matvec-krummhorn-2.data' using ($1+2):($4) title 'Intel Core 2 Duo' w l, \
     'matvec-datastar-2.data'     using ($1+2):($3) title 'IBM Power4'  w l, \
     'matvec-padoan-2.data'       using ($1+2):($2) title 'AMD Opteron' w l
