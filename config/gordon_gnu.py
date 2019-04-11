import os

is_arch_valid = 1

#python_lt_27 = 1

flags_arch = '-O3 -Wall -g'
#flags_arch = '-Wall -g'
flags_link  = '-rdynamic'

cc   = 'gcc'
f90  = 'gfortran'

#flags_prec_single = '-fdefault-real-4 -fdefault-double-8'
flags_prec_single = ''
flags_prec_double = '-fdefault-real-8 -fdefault-double-8'


libpath_fortran = ''
libs_fortran    = ['gfortran']

home = os.getenv('HOME')

hdf5_path    = os.environ['HDF5HOME']
hdf5_inc = hdf5_path + '/include'
hdf5_lib = hdf5_path + '/lib'

charm_path   = home + '/Charm/charm'
papi_inc = home + '/include'
papi_lib = home + '/lib'

png_path     = '/usr/lib64'
grackle_path = home + '/public/Grackle/src/clib'

serial_run   = charm_path + '/bin/charmrun  ++runscript ./run.sh +p1 '
parallel_run = charm_path + '/bin/charmrun  ++runscript ./run.sh +p4 '
parallel_arg = ''

