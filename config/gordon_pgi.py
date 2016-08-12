import os

is_arch_valid = 1

#python_lt_27 = 1

flags_arch = '-O3'
flags_link = '-pgf90libs'

cc   = 'pgcc'
f90  = 'pgf90'

libpath_fortran = ''
libs_fortran    = []

flags_prec_single = '-Mr4'
flags_prec_double = '-Mr8'

home = os.environ['HOME']

libpath_fortran = ''
libs_fortran    = ['gfortran']

home = os.environ["HOME"]
charm_path = home + '/Charm/charm'
papi_path  = home
hdf5_path  = os.environ['HDF5HOME']
hdf5_inc = hdf5_path + '/include'
hdf5_lib = hdf5_path + '/lib'
mpi_path   = os.environ['MPIHOME']

png_path   = '/usr/lib64'
grackle_path = home + '/Grackle/src/clib'
