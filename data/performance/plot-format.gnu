#========================================================================
set title  "Datastar (IBM Power4) Matvec Performance"
set xlabel "N (problem size = N^3)"
set ylabel "MFlop/s"

set output  "plot-format-datastar.png"
set terminal png

plot [][0:1800] \
     'datastar/matvec-O2.data' using ($1+2):($4) title 'Nonsymmetric' w l, \
     'datastar/matvec-O2.data' using ($1+2):($3) title 'Symmetric'  w l, \
     'datastar/matvec-O2.data' using ($1+2):($2) title 'Constant' w l

set output  "plot-format-datastar.eps"
set terminal postscript color eps font 24 size 6.5,4.5

plot [][0:1800] \
     'datastar/matvec-O2.data' using ($1+2):($4) title 'Nonsymmetric' w l, \
     'datastar/matvec-O2.data' using ($1+2):($3) title 'Symmetric'  w l, \
     'datastar/matvec-O2.data' using ($1+2):($2) title 'Constant' w l

#========================================================================
set title  "Kraken (AMD Opteron) Matvec Performance"
set xlabel "N (problem size = N^3)"
set ylabel "MFlop/s"

set output  "plot-format-kraken.png"
set terminal png

plot [][0:1800] \
     'kraken/matvec-O2.data' using ($1+2):($4) title 'Nonsymmetric' w l, \
     'kraken/matvec-O2.data' using ($1+2):($3) title 'Symmetric'  w l, \
     'kraken/matvec-O2.data' using ($1+2):($2) title 'Constant' w l

set output  "plot-format-kraken.eps"
set terminal postscript color eps font 24 size 6.5,4.5

plot [][0:1800] \
     'kraken/matvec-O2.data' using ($1+2):($4) title 'Nonsymmetric' w l, \
     'kraken/matvec-O2.data' using ($1+2):($3) title 'Symmetric'  w l, \
     'kraken/matvec-O2.data' using ($1+2):($2) title 'Constant' w l

#========================================================================
set title  "Padoan (AMD Opteron) Matvec Performance"
set xlabel "N (problem size = N^3)"
set ylabel "MFlop/s"

set output  "plot-format-padoan.png"
set terminal png

plot [][0:1800] \
     'padoan/matvec-O2.data' using ($1+2):($4) title 'Nonsymmetric' w l, \
     'padoan/matvec-O2.data' using ($1+2):($3) title 'Symmetric'  w l, \
     'padoan/matvec-O2.data' using ($1+2):($2) title 'Constant' w l

set output  "plot-format-padoan.eps"
set terminal postscript color eps font 24 size 6.5,4.5

plot [][0:1800] \
     'padoan/matvec-O2.data' using ($1+2):($4) title 'Nonsymmetric' w l, \
     'padoan/matvec-O2.data' using ($1+2):($3) title 'Symmetric'  w l, \
     'padoan/matvec-O2.data' using ($1+2):($2) title 'Constant' w l

#========================================================================
set title  "Intel Core 2 Q9450 Matvec Performance"
set xlabel "N (problem size = N^3)"
set ylabel "MFlop/s"

set output  "plot-format-intel-Q9450.png"
set terminal png

plot [][0:1800] \
     'intel-Q9450/matvec-O2.data' using ($1+2):($4) title 'Nonsymmetric' w l, \
     'intel-Q9450/matvec-O2.data' using ($1+2):($3) title 'Symmetric'  w l, \
     'intel-Q9450/matvec-O2.data' using ($1+2):($2) title 'Constant' w l

set output  "plot-format-intel-Q9450.eps"
set terminal postscript color eps font 24 size 6.5,4.5

plot [][0:1800] \
     'intel-Q9450/matvec-O2.data' using ($1+2):($4) title 'Nonsymmetric' w l, \
     'intel-Q9450/matvec-O2.data' using ($1+2):($3) title 'Symmetric'  w l, \
     'intel-Q9450/matvec-O2.data' using ($1+2):($2) title 'Constant' w l

#========================================================================


