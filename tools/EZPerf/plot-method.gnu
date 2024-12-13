set terminal png size 800,600
set termoption noenhanced
set output "plot-method.png"
set title "Enzo-E method times"
set xlabel "cycle"
set ylabel "time (s)"
set key bottom right
set grid
set log y
plot "method.data" using ($1):($2/1000000) title "total" w l lw 2 , \
     "method_comoving_expansion.data" using ($1):($2/1000000) title "comoving_expansion" w l lw 2 , \
     "method_gravity.data" using ($1):($2/1000000) title "gravity" w l lw 2 , \
     "method_pm_deposit.data" using ($1):($2/1000000) title "pm_deposit" w l lw 2 , \
     "method_pm_update.data" using ($1):($2/1000000) title "pm_update" w l lw 2 , \
     "method_ppm.data" using ($1):($2/1000000) title "ppm" w l lw 2 
