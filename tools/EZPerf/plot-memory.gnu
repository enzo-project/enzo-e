set terminal png size 800,600
set termoption noenhanced
set output "plot-memory.png"
set title "Enzo-E total memory"
set xlabel "cycle"
set ylabel "memory (GB)"
set key bottom right
set grid
set log y

plot "bytes-curr.data" using ($1):($2/1024/1024/1024) title "mem-curr" w l lw 2 , \
     "bytes-high.data" using ($1):($2/1024/1024/1024) title "mem-high" w l lw 2 , \
     "bytes-highest.data" using ($1):($2/1024/1024/1024) title "mem-highest" w l lw 2
