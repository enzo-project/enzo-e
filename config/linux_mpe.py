f90 = {}
cxx = {}
cc  = {}

is_arch_valid = 1

#flags_arch = '-g -Wall'
flags_arch = '-O3 -Wall'
flags_link = '-rdynamic -lpthread'

cc['charm']   = 'gcc'
cc['mpi']     = 'mpicc'
cc['serial']  = 'gcc'
cxx['mpi']    = 'mpic++'
cxx['serial'] = 'g++'
f90['charm']  = 'gfortran'
f90['mpi']    = 'gfortran'
f90['serial'] = 'gfortran'

libpath_fortran = ''
libs_fortran    = ['gfortran','lmpe','mpe']

charm_path  = '/home/bordner/charm/charm'
papi_path   = '/usr/local'
hdf5_path   = '/usr'
