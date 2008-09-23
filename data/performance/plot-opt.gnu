#========================================================================
set title  "Krummhorn Matvec Performance versus Optimization Level"
set xlabel "N (problem size = N^3)"
set ylabel "MFlop/s"

set output  "plot-opt-krummhorn.png"
set terminal png

plot [][0:1800] 'matvec-krummhorn-O0-2.data' using ($1+2):($4) title 'g++ -O0' w l, \
     'matvec-krummhorn-O1-2.data' using ($1+2):($4) title 'g++ -O1' w l, \
     'matvec-krummhorn-O2-2.data' using ($1+2):($4) title 'g++ -O2' w l, \
     'matvec-krummhorn-O3-2.data' using ($1+2):($4) title 'g++ -O3' w l

set output  "plot-opt-krummhorn.eps"
set terminal postscript color eps font 18

plot [][0:1800] 'matvec-krummhorn-O0-2.data' using ($1+2):($4) title 'g++ -O0' w l, \
     'matvec-krummhorn-O1-2.data' using ($1+2):($4) title 'g++ -O1' w l, \
     'matvec-krummhorn-O2-2.data' using ($1+2):($4) title 'g++ -O2' w l, \
     'matvec-krummhorn-O3-2.data' using ($1+2):($4) title 'g++ -O3' w l

#========================================================================
set title  "Padoan Matvec Performance versus Optimization Level"
set xlabel "N (problem size = N^3)"
set ylabel "MFlop/s"

set output  "plot-opt-padoan.png"
set terminal png

plot [][0:1800] 'matvec-padoan-O0-2.data' using ($1+2):($4) title 'g++ -O0' w l, \
     'matvec-padoan-O1-2.data' using ($1+2):($4) title 'g++ -O1' w l, \
     'matvec-padoan-O2-2.data' using ($1+2):($4) title 'g++ -O2' w l, \
     'matvec-padoan-O3-2.data' using ($1+2):($4) title 'g++ -O3' w l

set output  "plot-opt-padoan.eps"
set terminal postscript color eps font 18

plot [][0:1800] 'matvec-padoan-O0-2.data' using ($1+2):($4) title 'g++ -O0' w l, \
     'matvec-padoan-O1-2.data' using ($1+2):($4) title 'g++ -O1' w l, \
     'matvec-padoan-O2-2.data' using ($1+2):($4) title 'g++ -O2' w l, \
     'matvec-padoan-O3-2.data' using ($1+2):($4) title 'g++ -O3' w l

#========================================================================
set title  "Datastar Matvec Performance versus Optimization Level"
set xlabel "N (problem size = N^3)"
set ylabel "MFlop/s"

set output  "plot-opt-datastar.png"
set terminal png

plot [][0:1800] 'matvec-datastar-O0-2.data' using ($1+2):($4) title 'g++ -O0' w l, \
     'matvec-datastar-O2-2.data' using ($1+2):($4) title 'g++ -O2' w l, \
     'matvec-datastar-O3-2.data' using ($1+2):($4) title 'g++ -O3' w l, \
     'matvec-datastar-O4-2.data' using ($1+2):($4) title 'g++ -O4' w l, \
     'matvec-datastar-O5-2.data' using ($1+2):($4) title 'g++ -O5' w l

set output  "plot-opt-datastar.eps"
set terminal postscript color eps font 18

plot [][0:1800] 'matvec-datastar-O0-2.data' using ($1+2):($4) title 'g++ -O0' w l, \
     'matvec-datastar-O2-2.data' using ($1+2):($4) title 'g++ -O2' w l, \
     'matvec-datastar-O3-2.data' using ($1+2):($4) title 'g++ -O3' w l, \
     'matvec-datastar-O4-2.data' using ($1+2):($4) title 'g++ -O4' w l, \
     'matvec-datastar-O5-2.data' using ($1+2):($4) title 'g++ -O5' w l

#========================================================================
set title  "Stingmeyer Matvec Performance versus Optimization Level"
set xlabel "N (problem size = N^3)"
set ylabel "MFlop/s"

set output  "plot-opt-stingmeyer.png"
set terminal png

plot [][0:1800] 'matvec-stingmeyer-O0-2.data' using ($1+2):($4) title 'g++ -O0' w l, \
     'matvec-stingmeyer-O1-2.data' using ($1+2):($4) title 'g++ -O1' w l, \
     'matvec-stingmeyer-O2-2.data' using ($1+2):($4) title 'g++ -O2' w l, \
     'matvec-stingmeyer-O3-2.data' using ($1+2):($4) title 'g++ -O3' w l

set output  "plot-opt-stingmeyer.eps"
set terminal postscript color eps font 18

plot [][0:1800] 'matvec-stingmeyer-O0-2.data' using ($1+2):($4) title 'g++ -O0' w l, \
     'matvec-stingmeyer-O1-2.data' using ($1+2):($4) title 'g++ -O1' w l, \
     'matvec-stingmeyer-O2-2.data' using ($1+2):($4) title 'g++ -O2' w l, \
     'matvec-stingmeyer-O3-2.data' using ($1+2):($4) title 'g++ -O3' w l
#========================================================================
set title  "Kraken Matvec Performance versus Optimization Level"
set xlabel "N (problem size = N^3)"
set ylabel "MFlop/s"

set output  "plot-opt-kraken.png"
set terminal png

plot [][0:1800] 'matvec-kraken-O0-1.data' using ($1+2):($4) title 'g++ -O0' w l, \
     'matvec-kraken-O1-1.data' using ($1+2):($4) title 'g++ -O1' w l, \
     'matvec-kraken-O2-1.data' using ($1+2):($4) title 'g++ -O2' w l, \
     'matvec-kraken-O3-1.data' using ($1+2):($4) title 'g++ -O2' w l, \
     'matvec-kraken-O4-1.data' using ($1+2):($4) title 'g++ -O3' w l

set output  "plot-opt-kraken.eps"
set terminal postscript color eps font 18

plot [][0:1800] 'matvec-kraken-O0-1.data' using ($1+2):($4) title 'g++ -O0' w l, \
     'matvec-kraken-O1-1.data' using ($1+2):($4) title 'g++ -O1' w l, \
     'matvec-kraken-O2-1.data' using ($1+2):($4) title 'g++ -O2' w l, \
     'matvec-kraken-O3-1.data' using ($1+2):($4) title 'g++ -O2' w l, \
     'matvec-kraken-O4-1.data' using ($1+2):($4) title 'g++ -O3' w l
#========================================================================
