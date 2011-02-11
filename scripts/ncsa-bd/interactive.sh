#@ job_type = parallel
#@ node_usage = shared
#@ environment = COPY_ALL
#@ tasks_per_node = 4
#@ node = 1
#@ wall_clock_limit = 0:10:00
### uncomment below for a normal batch job
# #@ output = $(host).$(jobid).$(stepid).out
# #@ error = $(host).$(jobid).$(stepid).err

#@ queue
## uncomment for a normal batch job
# $HOME/a.out
