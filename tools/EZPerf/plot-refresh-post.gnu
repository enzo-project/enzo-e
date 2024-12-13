set terminal png size 800,600
set termoption noenhanced
set output "plot-refresh-post.png"
set title "Enzo-E refresh overhead"
set xlabel "cycle"
set ylabel "time (s)"
set key bottom right
set grid
set log y

plot "refresh_post.data" using ($1):($2/1000000) title "refresh total" w l lw 2 , \
     "refresh_recv_post.data" using ($1):($2/1000000) title "refresh_recv" w l lw 2 , \
     "refresh_child_post.data" using ($1):($2/1000000) title "refresh_child" w l lw 2 , \
     "refresh_exit_post.data" using ($1):($2/1000000) title "refresh_exit" w l lw 2
