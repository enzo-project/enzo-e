set title  "Matvec Performance versus Optimization Level"
set xlabel "N (problem size = N^3)"
set ylabel "MFlop/s"

set output  "plot-5.png"
set terminal png

plot 'matvec-krummhorn-O0-2.data' using ($1+2):($4) title 'g++ -O0' w l, \
     'matvec-krummhorn-O1-2.data' using ($1+2):($4) title 'g++ -O1' w l, \
     'matvec-krummhorn-O2-2.data' using ($1+2):($4) title 'g++ -O2' w l, \
     'matvec-krummhorn-O3-2.data' using ($1+2):($4) title 'g++ -O3' w l

set output  "plot-5.eps"
set terminal postscript color eps font 18

plot 'matvec-krummhorn-O0-2.data' using ($1+2):($4) title 'g++ -O0' w l, \
     'matvec-krummhorn-O1-2.data' using ($1+2):($4) title 'g++ -O1' w l, \
     'matvec-krummhorn-O2-2.data' using ($1+2):($4) title 'g++ -O2' w l, \
     'matvec-krummhorn-O3-2.data' using ($1+2):($4) title 'g++ -O3' w l
