# parses Enzo-p output and prints the following
#
#   number of solver iterations
#   time per iteration
#   number of blocks
#   number of particles
#   time per block
#   time per particle
#   time per zone
#   current memory
#   highest memory
#

BEGIN {
    block_zones_real = 1;
    block_zones_ghost = 1;
}
/ghost_depth/ { ghosts = $6;}
/iter 0000/  {time_start=$2;}
/final iter/ {
    max_iter=$7;
    printf ("%20s: %g\n", "solver max-iter",max_iter);
    time_stop = $2;
    time_per_iter = (time_stop-time_start)/max_iter; 
    printf ("%20s: %g\n", "time per iter",time_per_iter);
}
/Simulation cycle/ {cycle = $5;}
/bytes-high / {
    if (cycle == 1) {
	printf ("%20s: %g\n", "bytes-high",$6);
    }
}
/bytes-highest/ {
    if (cycle == 1) {
	printf ("%20s: %g\n", "bytes-highest",$6);
    }
}
/num-blocks/ {
    if (cycle == 1) {
	num_blocks = $6;
	printf ("%20s: %g\n", "num-blocks",num_blocks);
    }
}
/num-particles/ {
    if (cycle == 1) {
	num_particles = $6;
	printf ("%20s: %g\n","num-particles",num_particles);
    }
}
/root_blocks/ {
    n = $6;
    block_zones_real = block_zones_real * n;
    block_zones_ghost = block_zones_ghost * (2*ghosts * n);
}
/cycle 0002/ {
    printf ("%20s: %g\n", "time per block",time_per_iter/num_blocks);
    printf ("%20s: %g\n", "time per real zone",time_per_iter/(num_blocks*block_zones_real));
    printf ("%20s: %g\n", "time per zone",time_per_iter/(num_blocks*block_zones_ghost));
    printf ("%20s: %g\n", "time per particle",time_per_iter/num_particles);
}
/END CELLO/ {
    time_final = $2;
    printf ("%20s: %g\n", "Total time",time_final);
    printf ("%20s: %g\n", "Cycles",cycle);
    printf ("%20s: %g\n", "Avg time per cycle",(time_final-time_start)/cycle);
}

