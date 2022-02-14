import os

is_arch_valid = 1

#python_lt_27 = 1

flags_arch = '-O3 -Wall -g -ffast-math -funroll-loops' # -march=cascadelake
#flags_arch = '-Wall -g'
flags_link  = '-rdynamic'

cc   = 'gcc'
f90  = 'gfortran'

#flags_prec_single = '-fdefault-real-4 -fdefault-double-8'
flags_prec_single = ''
flags_prec_double = '-fdefault-real-8 -fdefault-double-8'

libpath_fortran = ''
libs_fortran    = ['gfortran']
flags_arch_fortran = '-ffixed-line-length-132'

home = os.environ['HOME']

hdf5_path    = os.environ['HDF5_HOME']
hdf5_inc = hdf5_path + '/include'
hdf5_lib = hdf5_path + '/lib'

boost_path = os.environ['BOOST_HOME']
boost_inc = boost_path 
boost_lib = boost_path 

charm_path = os.environ['CHARM_HOME']

papi_inc = home + '/include'
papi_lib = home + '/lib'

png_path     = '/usr/lib64'
grackle_path = home + '/local'

cello_var = os.environ.get('CELLO_VAR',"mpi-gcc")

if (cello_var == "mpi-gcc"):
    parallel_run = charm_path + "/bin/charmrun +p4 "
    parallel_arg = " "
    smp = 0
elif (cello_var == "mpi-gcc-smp"):
    parallel_run = charm_path + "/bin/charmrun +p4 "
    parallel_arg = " ++processPerHost 1 ++ppn 2 "
    smp = 1
elif (cello_var == "net-gcc"):
    parallel_run = charm_path + "/bin/charmrun +p4 ++local "
    parallel_arg = " "
    smp = 0
elif (cello_var == "net-gcc-smp"):
    parallel_run = charm_path + "/bin/charmrun +p4 ++local "
    parallel_arg = " ++processPerHost 1 ++ppn 2 "
    smp = 1
elif (cello_var == "ucx-gcc"):
    parallel_run = charm_path + "/bin/charmrun +p4 "
    parallel_arg = " "
    smp = 0
elif (cello_var == "ucx-gcc-smp"):
    parallel_run = charm_path + "/bin/charmrun +p4 "
    parallel_arg = " ++processPerHost 1 ++ppn 2 "
    smp = 1



