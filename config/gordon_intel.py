import os

is_arch_valid = 1

#flags_arch = '-g'
#flags_arch = '-Ofast' # ERROR: ipo linking
flags_arch = '-O3'
flags_link  = ''

cc   = 'icc'
f90  = 'ifort'

flags_prec_single = '-real-size-4 -double-size-8'
flags_prec_double = '-real-size-8 -double-size-8'

libpath_fortran = ''
libs_fortran    = ['imf','ifcore','ifport','stdc++','intlc']

home = os.environ['HOME']

libpath_fortran = ''
libs_fortran    = ['ifcore']

home = os.environ["HOME"]
charm_path = home + '/Charm/charm'
papi_path  = home
hdf5_path  = os.environ['HDF5HOME']
mpi_path   = os.environ['MPIHOME']

png_path   = '/usr/lib64'
grackle_path = home + '/Grackle/src/clib'

