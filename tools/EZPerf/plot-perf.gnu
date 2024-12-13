set terminal png size 800,600
set termoption noenhanced
set output "plot-perf.png"
set title "Enzo-E performance summary"
set xlabel "cycle"
set ylabel "time (s)"
set key bottom right
set grid
set log y
plot 'cycle.data' using ($1):($2) title 'cycle' w l lw 2 , \
     'method.data' using ($1):($2/1000000) title 'method' w l lw 2 , \
     'solver.data' using ($1):($2/1000000) title 'solver' w l lw 2 , \
     'refresh.data' using ($1):($2/1000000) title 'refresh' w l lw 2, \
     'adapt.data' using ($1):($2/1000000) title 'adapt' w l lw 2 , \
     'adapt_post.data' using ($1):($2/1000000) title 'adapt_post' w l lw 2

