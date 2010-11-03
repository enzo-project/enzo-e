#========================================================================
set title  "Datastar (IBM Power4) Matvec Performance"
set xlabel "N (problem size = N^3)"
set ylabel "MFlop/s"

set output  "plot-format-datastar.png"
set terminal png

plot [][0:1800] 'matvec-datastar-O2-2.data' using ($1+2):($4) title 'Nonsymmetric' w l, \
            'matvec-datastar-O2-2.data' using ($1+2):($3) title 'Symmetric'  w l, \
            'matvec-datastar-O2-2.data' using ($1+2):($2) title 'Constant' w l

set output  "plot-format-datastar.eps"
set terminal postscript color eps font 24 size 6.5,4.5

plot [][0:1800] 'matvec-datastar-O2-2.data' using ($1+2):($4) title 'Nonsymmetric' w l, \
            'matvec-datastar-O2-2.data' using ($1+2):($3) title 'Symmetric'  w l, \
            'matvec-datastar-O2-2.data' using ($1+2):($2) title 'Constant' w l

#========================================================================
set title  "Kraken (AMD Opteron) Matvec Performance"
set xlabel "N (problem size = N^3)"
set ylabel "MFlop/s"

set output  "plot-format-kraken.png"
set terminal png

plot [][0:1800] 'matvec-kraken-O2-1.data' using ($1+2):($4) title 'Nonsymmetric' w l, \
            'matvec-kraken-O2-1.data' using ($1+2):($3) title 'Symmetric'  w l, \
            'matvec-kraken-O2-1.data' using ($1+2):($2) title 'Constant' w l

set output  "plot-format-kraken.eps"
set terminal postscript color eps font 24 size 6.5,4.5

plot [][0:1800] 'matvec-kraken-O2-1.data' using ($1+2):($4) title 'Nonsymmetric' w l, \
            'matvec-kraken-O2-1.data' using ($1+2):($3) title 'Symmetric'  w l, \
            'matvec-kraken-O2-1.data' using ($1+2):($2) title 'Constant' w l

#========================================================================
set title  "Padoan (AMD Opteron) Matvec Performance"
set xlabel "N (problem size = N^3)"
set ylabel "MFlop/s"

set output  "plot-format-padoan.png"
set terminal png

plot [][0:1800] 'matvec-padoan-O2-2.data' using ($1+2):($4) title 'Nonsymmetric' w l, \
            'matvec-padoan-O2-2.data' using ($1+2):($3) title 'Symmetric'  w l, \
            'matvec-padoan-O2-2.data' using ($1+2):($2) title 'Constant' w l

set output  "plot-format-padoan.eps"
set terminal postscript color eps font 24 size 6.5,4.5

plot [][0:1800] 'matvec-padoan-O2-2.data' using ($1+2):($4) title 'Nonsymmetric' w l, \
            'matvec-padoan-O2-2.data' using ($1+2):($3) title 'Symmetric'  w l, \
            'matvec-padoan-O2-2.data' using ($1+2):($2) title 'Constant' w l

#========================================================================
set title  "Krummhorn (Intel Core 2 Duo) Matvec Performance"
set xlabel "N (problem size = N^3)"
set ylabel "MFlop/s"

set output  "plot-format-krummhorn.png"
set terminal png

plot [][0:1800] 'matvec-krummhorn-O2-2.data' using ($1+2):($4) title 'Nonsymmetric' w l, \
            'matvec-krummhorn-O2-2.data' using ($1+2):($3) title 'Symmetric'  w l, \
            'matvec-krummhorn-O2-2.data' using ($1+2):($2) title 'Constant' w l

set output  "plot-format-krummhorn.eps"
set terminal postscript color eps font 24 size 6.5,4.5

plot [][0:1800] 'matvec-krummhorn-O2-2.data' using ($1+2):($4) title 'Nonsymmetric' w l, \
            'matvec-krummhorn-O2-2.data' using ($1+2):($3) title 'Symmetric'  w l, \
            'matvec-krummhorn-O2-2.data' using ($1+2):($2) title 'Constant' w l

#========================================================================
set title  "Stingmeyer (Intel Core 2 Quad) Matvec Performance"
set xlabel "N (problem size = N^3)"
set ylabel "MFlop/s"

set output  "plot-format-stingmeyer.png"
set terminal png

plot [][0:1800] 'matvec-stingmeyer-O2-2.data' using ($1+2):($4) title 'Nonsymmetric' w l, \
            'matvec-stingmeyer-O2-2.data' using ($1+2):($3) title 'Symmetric'  w l, \
            'matvec-stingmeyer-O2-2.data' using ($1+2):($2) title 'Constant' w l

set output  "plot-format-stingmeyer.eps"
set terminal postscript color eps font 24 size 6.5,4.5

plot [][0:1800] 'matvec-stingmeyer-O2-2.data' using ($1+2):($4) title 'Nonsymmetric' w l, \
            'matvec-stingmeyer-O2-2.data' using ($1+2):($3) title 'Symmetric'  w l, \
            'matvec-stingmeyer-O2-2.data' using ($1+2):($2) title 'Constant' w l

#========================================================================


