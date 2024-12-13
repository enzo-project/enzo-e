set terminal png size 800,600
set termoption noenhanced
set output "plot-mesh.png"
set title "Enzo-E blocks per level"
set xlabel "cycle"
set ylabel "blocks"
set key bottom right
set grid
set log y

plot "total-blocks.data" using ($1):($2) title "total-blocks" w l lw 2 , \
     "leaf-blocks.data" using ($1):($2) title "leaf-blocks" w l lw 2 , \
     "blocks-level_0.data" using ($1):($2) title "blocks-level_0" w l lw 2 , \
     "blocks-level_1.data" using ($1):($2) title "blocks-level_1" w l lw 2 , \
     "blocks-level_2.data" using ($1):($2) title "blocks-level_2" w l lw 2 , \
     "blocks-level_3.data" using ($1):($2) title "blocks-level_3" w l lw 2 , \
     "blocks-level_4.data" using ($1):($2) title "blocks-level_4" w l lw 2,
     "blocks-level_5.data" using ($1):($2) title "blocks-level_5" w l lw 2,
     "blocks-level_6.data" using ($1):($2) title "blocks-level_6" w l lw 2
