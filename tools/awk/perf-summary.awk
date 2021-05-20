#!/usr/bin/awk -f

# parses Enzo-E output and prints the following
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
    max_level = 0;
    min_level = 0;
    done = 0;
    root_blocks = 1;
    root_size = 1;
    total_root_size = 1;
    i_root_blocks = 0;
    i_root_size = 0;
    bytes_high_start = 0;
    num_blocks_start = 0;
    zones_per_block = 1;
    total_zones_per_block = 1;
    first_solver = 1;
    index_solver = 0;
    num_solvers = 0;
}

/DEBUG/ { next }

# Solver cg iteration rate
/Solver / {
    time = $2;
    name = $4;
    iters = $6;

    if (iters == 0) {
	time_cycle_start[name]=time;
	iters_prev[name] = 0;
    } else {
	time_total[name]=time_total[name]+(time-time_cycle_start[name]);
	iters_cycle[name] = iters - iters_prev[name]
	iters_prev[name] = iters;
	iters_total[name] = iters_total[name]+iters_cycle[name];
	time_cycle_start[name]=time;
    }
}
/processors/{num_processors = $6; }
# Block size

/ghost_depth/ { ghosts = $6;}

/root_blocks/ {
    if ($7 ~ /^[0-9]+$/) {
	root_blocks3[0] = $7
	root_blocks = root_blocks * $7;
    } else root_blocks3[0] = 1;
    if ($8 ~ /^[0-9]+$/) {
	root_blocks3[1] = $8
	root_blocks = root_blocks * $8;
    } else root_blocks3[1] = 1;
    if ($9 ~ /^[0-9]+$/) {
	root_blocks3[2] = $9
	root_blocks = root_blocks * $9;
    } else root_blocks3[2] = 1;
}

/root_size/ {
    if ($7 ~ /^[0-9]+$/) {
	root_size = root_size * $7;
	total_root_size = total_root_size * ($7+2*ghosts)
    }
    if ($8 ~ /^[0-9]+$/) {
	root_size = root_size * $8;
	total_root_size = total_root_size * ($8+2*ghosts)
    }
    if ($9 ~ /^[0-9]+$/) {
	root_size = root_size * $9;
	total_root_size = total_root_size * ($9+2*ghosts)
    }
}


/Simulation cycle 0000/ {
    cycle_start=$2;
}

/Simulation cycle / {
    cycle_last=$2;
    cycle = $5;
}

# Linear Solver
# / iter /  {
#     # ASSUMES ONLY FIRST AND LAST ITERATIONS OUTPUT
#     # (i.e. monitor_iter = 0 for all solvers)
#     time=$2
#     name=$4;
#     iter=$6;
#     if (iter == "0000") {
# 	index_solver++;
# 	if (index_solver > num_solvers) num_solvers = index_solver;
#     } else {
# 	index_solver--;
#     }
#     solver_time[name] = time - solver_time[name];
# }

# Linear Solver
/final iter/ {
    max_iter      = $7;
    time_stop     = $2;
    time_per_iter = (time_stop-time_start)/max_iter; 
}

/bytes-high / {
    bytes_high = $6;
    if (bytes_high_start == 0) bytes_high_start = $6;
}

/bytes-highest/ {  bytes_highest = $6; }

/num-leaf-blocks/ {
    if (num_blocks_start == 0) num_blocks_start = $6;
    if ($6 > num_blocks) num_blocks = $6;
    num_blocks_total = num_blocks_total + $6;
}

/num-particles/ {
    num_particles = $7;
    if (num_particles_start == 0) {
	num_particles_start = num_particles;
    }
}
/Mesh:root_size/ {
    if ($7 ~ /^[0-9]+$/) {
	root_size3[0] = $7;
    }
    if ($8 ~ /^[0-9]+$/) {
	root_size3[1] = $8;
    }
    if ($9 ~ /^[0-9]+$/) {
	root_size3[2] = $9;
    }
}
/Adapt:max_level/ { max_level = $9 }
/Adapt:min_level/ { min_level = $9 }

/Mesh:root_rank/ {root_rank = $6; }

/END ENZO/ {
    done = 1;
    time_final = $2;
}

/simulation time-usec/  {time_simulation  = $6}
/simulation PAPI_FP_INS/ {
    if (fp_ins_simulation == 0) {
	fp_ins_simulation_start = $6;
    }
    fp_ins_simulation = $6;
}
/simulation PAPI_FP_OPS/ {
    if (fp_ops_simulation == 0) {
	fp_ops_simulation_start = $6;
    }
    fp_ops_simulation = $6;
}
/initial time-usec/  {time_initial  = $6}
/cycle time-usec/    {time_cycle    = $6}
/compute time-usec/  {time_compute  = $6};
/adapt_apply time-usec/    {time_adapt_apply    = $6};
/adapt_apply_sync time-usec/    {time_adapt_apply_sync    = $6};
/adapt_notify time-usec/    {time_adapt_notify    = $6};
/adapt_notify_sync time-usec/    {time_adapt_notify_sync    = $6};
/adapt_update time-usec/    {time_adapt_update    = $6};
/adapt_update_sync time-usec/    {time_adapt_update_sync    = $6};
/adapt_end time-usec/    {time_adapt_end    = $6};
/adapt_end_sync time-usec/    {time_adapt_end_sync    = $6};
/refresh_store time-usec/  {time_refresh_store  = $6};
/refresh_child time-usec/  {time_refresh_child  = $6};
/refresh_exit time-usec/  {time_refresh_exit  = $6};
/refresh_store_sync time-usec/  {time_refresh_store_sync  = $6};
/refresh_child_sync time-usec/  {time_refresh_child_sync  = $6};
/refresh_exit_sync time-usec/  {time_refresh_exit_sync  = $6};
/output time-usec/   {time_output   = $6};
/stopping time-usec/ {time_stopping = $6};

END {

    zones_per_block       = root_size / root_blocks;
    total_zones_per_block = total_root_size  / root_blocks;
    if (cycle != 0) {
	time_per_cycle    = (cycle_last-cycle_start)/cycle;
    }
    if (time_per_cycle != 0) {
	zone_rate         = zones_per_block/time_per_cycle*num_blocks;
    }
    if (zone_rate != 0) {
	zone_us           = 1.0e6 / zone_rate*num_processors ;
    }
    
    format = "%25s: %8.4f\n";
    format_ind = "   %25s: %8.4f\n";
    format_int = "%25s: %8ld\n";
    format2 = "%25s: %8.4f [%8.4f ]\n";
    format2_ind = "   %25s: %8.4f [%8.4f ]\n";
    format3 = "%25s: %4d %d %d\n";
    format0 = "%25s\n";

    if (done == 1) {
	print ("Run completed");
    } else {
	print ("RUN INCOMPLETE");
    }

    print "\nSUMMARY\n"
    
    printf (format, "Total time (s)",time_final);
    printf (format_int, "Processes", num_processors);
    printf (format_int, "Cycles",cycle);
    printf (format, "Init time (s)",cycle_start);
    printf (format, "Zone rate (Mz/s)",1.0e-6*zone_rate);
    printf (format, "Zone update time (µs/z/p)",zone_us);
    printf (format, "Time per cycle (s)",time_per_cycle);
    if (fp_ins_simulation > 0) {
	printf (format_int, "GFlops", fp_ins_simulation);
	printf (format, "GFlops per cycle", (fp_ins_simulation-fp_ins_simulation_start)/cycle);
	printf (format, "GFlop ins. rate", 1e-9*fp_ins_simulation/time_final);
    }
    if (fp_ops_simulation > 0) {
	printf (format_int, "GFlops", fp_ops_simulation);
	printf (format, "GFlops per cycle", (fp_ops_simulation-fp_ops_simulation_start)/cycle);
	printf (format, "GFlop ops. rate", 1e-9*fp_ops_simulation/time_final);
    }

    printf "\nHIERARCHY\n"
    
    printf (format_int, "rank",root_rank);
    printf (format3,"root blocks",root_blocks3[0],root_blocks3[1],root_blocks3[2]);
    printf (format3,"root size",root_size3[0],root_size3[1],root_size3[2]);
    printf (format3,"effective size",
	    root_size3[0]*2^max_level,
	    root_size3[1]*2^max_level,
	    root_size3[2]*2^max_level);
    bx=(root_size3[0]/root_blocks3[0]);
    by=(root_size3[1]/root_blocks3[1]);
    bz=(root_size3[2]/root_blocks3[2]);
    printf (format3,"block size",bx,by,bz);
    printf (format_int, "max-level",max_level);
    printf (format_int, "min-level",min_level);
    printf (format_int, "num-blocks",num_blocks);
    printf (format_int, "num-zones",num_blocks*bx*by*bz);
#    printf (format, "blocks change",num_blocks / num_blocks_start);

    printf ("\nMEMORY\n");
    
    scale = 1;
    prefix="";

    if (bytes_high > 1e3) {
	scale = scale * 1024;
	prefix="K";
    }
    if (bytes_high > 1e6) {
	scale = scale * 1024;
	prefix="M";
    }
    if (bytes_high > 1e9) {
	scale = scale * 1024;
	prefix="G";
    }
    if (bytes_high > 1e12) {
	scale = scale * 1024;
	prefix="T";
    }

    printf (format, "bytes-high (" prefix "B)",bytes_high/scale);
    printf (format, "bytes-highest (" prefix "B)", bytes_highest/scale);
    if (bytes_high_start != 0) {
	printf (format, "memory change",bytes_high / bytes_high_start);
    }
    printf (format, "bytes per real zone",bytes_high / (num_blocks*zones_per_block));
    printf (format, "bytes per zone",bytes_high / (num_blocks*total_zones_per_block));


    printf "\nPHASES\n"

    t_scale=1.0/num_processors*0.000001;
    printf (format, "Simulation",time_simulation*t_scale);
    printf (format, "Initial",time_initial*t_scale);
    printf (format, "Cycling",time_cycle*t_scale);
    printf("\n");
    printf (format, "Compute",time_compute*t_scale);
    time_adapt = time_adapt_apply + time_adapt_update + time_adapt_notify + time_adapt_end;
    time_adapt_sync = time_adapt_apply_sync + time_adapt_update_sync + time_adapt_notify_sync + time_adapt_end_sync;
    printf (format2, "Adapt",
	    t_scale*time_adapt,
	    t_scale*time_adapt_sync);
    printf (format2_ind, "apply",
	    t_scale*time_adapt_apply,
	    t_scale*time_adapt_apply_sync);
    printf (format2_ind, "update",
	    t_scale*time_adapt_update,
	    t_scale*time_adapt_update_sync);
    printf (format2_ind, "notify",
	    t_scale*time_adapt_notify,
	    t_scale*time_adapt_notify_sync);
    printf (format2_ind, "end",
	    t_scale*time_adapt_end,
	    t_scale*time_adapt_end_sync);
    time_refresh = time_refresh_store + time_refresh_child + time_refresh_exit;
    time_refresh_sync = time_refresh_store_sync + time_refresh_child_sync + time_refresh_exit_sync;
    printf (format2, "Refresh",
	    t_scale*time_refresh,
	    t_scale*time_refresh_sync);

    printf (format2_ind, "store",
	    t_scale*time_refresh_store,
	    t_scale*time_refresh_store_sync);
    printf (format2_ind, "child",
	    t_scale*time_refresh_child,
	    t_scale*time_refresh_child_sync);
    printf (format2_ind, "exit",
	    t_scale*time_refresh_exit,
	    t_scale*time_refresh_exit_sync);
    
    printf (format, "Output",time_output*t_scale);
    printf (format, "Stopping",time_stopping*t_scale);
    

    printf ("\nSOLVER\n");

    printf ("\n   Solver time\n");
    for (name in time_total) {
	printf (format_ind, "solver-time-"name, time_total[name]);
    }
    printf ("\n   Solver rate (i/b/p)\n");
    for (name in time_total) {
	rate_avg[name] = iters_total[name]/time_total[name]*num_blocks_total/num_processors;

	printf (format_ind, "solver-rate-" name,rate_avg[name]);
    }
	
    
#    if (time_per_iter > 0) {
#	print;
#	printf (format, "time per block",time_per_iter/num_blocks);
#	printf (format, "time per real zone",time_per_iter/(num_blocks*zones_per_block));
#	printf (format, "time per zone",time_per_iter/(num_blocks*total_zones_per_block));
#	printf (format_int, "solver max-iter",max_iter);
#	printf (format, "time per iter",time_per_iter);
#    } else {
#	print " not called"
#    }

    printf ("\nPARTICLES");
    
    if (num_particles != 0) {
	printf ("\n");
	printf (format_int,"num-particles",num_particles);
	printf (format, "time per particle",time_per_iter/num_particles);
	printf (format, "particles change",num_particles / num_particles_start);
    } else {
	print " not used"
    }
}

