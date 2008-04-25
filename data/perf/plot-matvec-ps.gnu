set output  "matvec.eps"
set terminal postscript color eps

plot 'matvec-datastar-1.data'  using ($1+2):($3) title 'IBM Power4' w l, \
     'matvec-padoan-1.data'    using ($1+2):($3) title 'AMD Opteron' w l, \
     'matvec-krummhorn-1.data' using ($1+2):($3) title 'Intel Core2 Duo' w l
