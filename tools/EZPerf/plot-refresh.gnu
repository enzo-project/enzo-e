set terminal png size 800,600
set termoption noenhanced
set output "plot-refresh.png"
set title "Enzo-E refresh time"
set xlabel "cycle"
set ylabel "time (s)"
set key bottom right
set grid
set log y

plot "refresh.data" using ($1):($2/1000000) title "refresh total" w l lw 2 , \
     "refresh_recv.data" using ($1):($2/1000000) title "refresh_recv" w l lw 2 , \
     "refresh_child.data" using ($1):($2/1000000) title "refresh_child" w l lw 2 , \
     "refresh_exit.data" using ($1):($2/1000000) title "refresh_exit" w l lw 2
