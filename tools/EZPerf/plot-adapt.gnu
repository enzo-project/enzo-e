set terminal png size 800,600
set termoption noenhanced
set output "plot-adapt.png"
set title "Enzo-E adapt time"
set xlabel "cycle"
set ylabel "time (s)"
set key bottom right
set grid
set log y

plot "adapt.data" using ($1):($2/1000000) title "adapt total" w l lw 2 , \
     "adapt_called.data" using ($1):($2/1000000) title "adapt_called" w l lw 2 , \
     "adapt_enter.data" using ($1):($2/1000000) title "adapt_enter" w l lw 2 , \
     "adapt_exit.data" using ($1):($2/1000000) title "adapt_exit" w l lw 2 , \
     "adapt_next.data" using ($1):($2/1000000) title "adapt_next" w l lw 2 , \
     "adapt_recv_level.data" using ($1):($2/1000000) title "adapt_recv_level" w l lw 2 , \
     "adapt_update.data" using ($1):($2/1000000) title "adapt_update" w l lw 2 , \
