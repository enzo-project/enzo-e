import os

is_arch_valid = 1

flags_arch = '-O3 -Wall -g'
flags_link  = '-rdynamic'

cc   = 'gcc'
f90  = 'gfortran'

flags_prec_single = '-fdefault-real-4 -fdefault-double-8'
flags_prec_double = '-fdefault-real-8 -fdefault-double-8'

home = os.environ['HOME']

libpath_fortran = ''
libs_fortran    = ['gfortran']

home = os.environ["HOME"]
charm_path = home + '/Charm/charm'
papi_path  = home
#hdf5_path  = os.environ['HDF5HOME']
hdf5_path  = home + '/public'
#mpi_path   = os.environ['MPIHOME']
mpi_path   = home + '/public'

png_path   = '/usr/lib64'
grackle_path = home + '/Grackle/src/clib'
