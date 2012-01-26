f90 = {}
cxx = {}
cc  = {}

is_arch_valid = 1

flags_arch = '-g -Wall -O'
flags_link = '-rdynamic'

cc['mpi']     = 'mpicc'
cc['serial']  = 'gcc'
cxx['mpi']    = 'mpic++'
cxx['serial'] = 'g++'
f90['charm']  = 'gfortran'
f90['mpi']    = 'gfortran'
f90['serial'] = 'gfortran'

libpath_fortran = ''
libs_fortran    = ['gfortran']

charm_path  = '/home/bordner/charm/charm'
papi_path   = '/usr/local'
hdf5_path   = '/usr'
