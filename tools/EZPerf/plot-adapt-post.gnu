set terminal png size 800,600
set termoption noenhanced
set output "plot-adapt-post.png"
set title "Enzo-E adapt overhead"
set xlabel "cycle"
set ylabel "time (s)"
set key bottom right
set grid
set log y

plot "adapt_post.data" using ($1):($2/1000000) title "adapt total" w l lw 2 , \
     "adapt_called_post.data" using ($1):($2/1000000) title "adapt_called" w l lw 2 , \
     "adapt_enter_post.data" using ($1):($2/1000000) title "adapt_enter" w l lw 2 , \
     "adapt_exit_post.data" using ($1):($2/1000000) title "adapt_exit" w l lw 2 , \
     "adapt_next_post.data" using ($1):($2/1000000) title "adapt_next" w l lw 2 , \
     "adapt_recv_level_post.data" using ($1):($2/1000000) title "adapt_recv_level" w l lw 2 , \
     "adapt_update_post.data" using ($1):($2/1000000) title "adapt_update" w l lw 2 , \
