set terminal png size 800,600
set termoption noenhanced
set output "plot-balance.png"
set title "Enzo-E load balance efficiency"
set xlabel "cycle"
set ylabel "efficiency (avg-time/max-time)"
set key bottom right
set grid

plot "eff-blocks-core.data" using ($1):($2) title "zones cores" w l lw 2 , \
     "eff-blocks-node.data" using ($1):($2) title "zones nodes" w l lw 2 , \
     "eff-particles-core.data" using ($1):($2) title "particles cores" w l lw 2 , \
     "eff-particles-node.data" using ($1):($2) title "particles nodes" w l lw 2 , \

