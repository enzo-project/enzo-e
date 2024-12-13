set terminal png size 800,600
set termoption noenhanced
set output "plot-solver.png"
set title "Enzo-E solver times"
set xlabel "cycle"
set ylabel "time (s)"
set key bottom right
set grid
set log y
plot 'solver.data' using ($1):($2/1000000) title 'total' w l lw 2 , \
     'solver_dd.data' using ($1):($2/1000000) title 'dd' w l lw 2 , \
     'solver_dd_domain.data' using ($1):($2/1000000) title 'dd_domain' w l lw 2 , \
     'solver_dd_root.data' using ($1):($2/1000000) title 'dd_root' w l lw 2 , \
     'solver_dd_smooth.data' using ($1):($2/1000000) title 'dd_smooth' w l lw 2 , \
     'solver_root_coarse.data' using ($1):($2/1000000) title 'root_coarse' w l lw 2 , \
     'solver_root_post.data' using ($1):($2/1000000) title 'root_post' w l lw 2 , \
     'solver_root_pre.data' using ($1):($2/1000000) title 'root_pre' w l lw 2 

