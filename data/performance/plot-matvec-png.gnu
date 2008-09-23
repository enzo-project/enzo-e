set output  "matvec.png"
set terminal png

plot 'matvec-datastar-1.data'  using ($1+2):($3) title 'IBM Power4' w l, \
     'matvec-padoan-1.data'    using ($1+2):($3) title 'AMD Opteron' w l, \
     'matvec-krummhorn-1.data' using ($1+2):($3) title 'Intel Core2 Duo 32-bit OS' w l, \
     'matvec-krummhorn-64-1.data' using ($1+2):($3) title 'Intel Core2 Duo 64-bit OS' w l

