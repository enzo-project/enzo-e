#========================================================================
set title  "Datastar (IBM Power4) Matvec Performance"
set xlabel "N (problem size = N^3)"
set ylabel "MFlop/s"

set output  "plot-opt-datastar.png"
set terminal png

plot [][0:1800] \
     'datastar/matvec-O0.data' using ($1+2):($4) title 'g++ -O0' w l, \
     'datastar/matvec-O2.data' using ($1+2):($4) title 'g++ -O2' w l, \
     'datastar/matvec-O3.data' using ($1+2):($4) title 'g++ -O3' w l, \
     'datastar/matvec-O4.data' using ($1+2):($4) title 'g++ -O4' w l, \
     'datastar/matvec-O5.data' using ($1+2):($4) title 'g++ -O5' w l

set output  "plot-opt-datastar.eps"
set terminal postscript color eps font 24 size 6.5,4.5

plot [][0:1800] \
     'datastar/matvec-O0.data' using ($1+2):($4) title 'g++ -O0' w l, \
     'datastar/matvec-O2.data' using ($1+2):($4) title 'g++ -O2' w l, \
     'datastar/matvec-O3.data' using ($1+2):($4) title 'g++ -O3' w l, \
     'datastar/matvec-O4.data' using ($1+2):($4) title 'g++ -O4' w l, \
     'datastar/matvec-O5.data' using ($1+2):($4) title 'g++ -O5' w l

#========================================================================
set title  "Kraken (AMD Opteron) Matvec Performance"
set xlabel "N (problem size = N^3)"
set ylabel "MFlop/s"

set output  "plot-opt-kraken.png"
set terminal png

plot [][0:1800] \
     'kraken/matvec-O0.data' using ($1+2):($4) title 'pgCC -O0' w l, \
     'kraken/matvec-O1.data' using ($1+2):($4) title 'pgCC -O1' w l, \
     'kraken/matvec-O2.data' using ($1+2):($4) title 'pgCC -O2' w l, \
     'kraken/matvec-O3.data' using ($1+2):($4) title 'pgCC -O3' w l, \
     'kraken/matvec-O4.data' using ($1+2):($4) title 'pgCC -O4' w l

set output  "plot-opt-kraken.eps"
set terminal postscript color eps font 24 size 6.5,4.5 

plot [][0:1800] \
     'kraken/matvec-O0.data' using ($1+2):($4) title 'pgCC -O0' w l, \
     'kraken/matvec-O1.data' using ($1+2):($4) title 'pgCC -O1' w l, \
     'kraken/matvec-O2.data' using ($1+2):($4) title 'pgCC -O2' w l, \
     'kraken/matvec-O3.data' using ($1+2):($4) title 'pgCC -O3' w l, \
     'kraken/matvec-O4.data' using ($1+2):($4) title 'pgCC -O4' w l
#========================================================================
set title  "Padoan (AMD Opteron) Matvec Performance"
set xlabel "N (problem size = N^3)"
set ylabel "MFlop/s"

set output  "plot-opt-padoan.png"
set terminal png

plot [][0:1800] \
     'padoan/matvec-O0.data' using ($1+2):($4) title 'g++ -O0' w l, \
     'padoan/matvec-O1.data' using ($1+2):($4) title 'g++ -O1' w l, \
     'padoan/matvec-O2.data' using ($1+2):($4) title 'g++ -O2' w l, \
     'padoan/matvec-O3.data' using ($1+2):($4) title 'g++ -O3' w l

set output  "plot-opt-padoan.eps"
set terminal postscript color eps font 24 size 6.5,4.5 

plot [][0:1800] \
     'padoan/matvec-O0.data' using ($1+2):($4) title 'g++ -O0' w l, \
     'padoan/matvec-O1.data' using ($1+2):($4) title 'g++ -O1' w l, \
     'padoan/matvec-O2.data' using ($1+2):($4) title 'g++ -O2' w l, \
     'padoan/matvec-O3.data' using ($1+2):($4) title 'g++ -O3' w l

#========================================================================
set title  "Intel Core 2 Q9450 Matvec Performance"
set xlabel "N (problem size = N^3)"
set ylabel "MFlop/s"

set output  "plot-opt-intel-Q9450.png"
set terminal png

plot [][0:1800] \
     'intel-Q9450/matvec-O0.data' using ($1+2):($4) title 'g++ -O0' w l, \
     'intel-Q9450/matvec-O1.data' using ($1+2):($4) title 'g++ -O1' w l, \
     'intel-Q9450/matvec-O2.data' using ($1+2):($4) title 'g++ -O2' w l, \
     'intel-Q9450/matvec-O3.data' using ($1+2):($4) title 'g++ -O3' w l

set output  "plot-opt-intel-Q9450.eps"
set terminal postscript color eps font 24 size 6.5,4.5 

plot [][0:1800] \
     'intel-Q9450/matvec-O0.data' using ($1+2):($4) title 'g++ -O0' w l, \
     'intel-Q9450/matvec-O1.data' using ($1+2):($4) title 'g++ -O1' w l, \
     'intel-Q9450/matvec-O2.data' using ($1+2):($4) title 'g++ -O2' w l, \
     'intel-Q9450/matvec-O3.data' using ($1+2):($4) title 'g++ -O3' w l
#========================================================================
