f90 = {}
cxx = {}
cc  = {}

is_arch_valid = 1

# flags_arch = '-g -Wall'
flags_arch = '-O3'
flags_link  = '-rdynamic'

cc['mpi']     = 'mpicc'
cc['serial']  = 'gcc'
cxx['mpi']    = 'mpiCC'
cxx['serial'] = 'g++'
f90['charm']  = 'gfortran'
f90['mpi']    = 'gfortran'
f90['serial'] = 'gfortran'

libpath_fortran = ''
libs_fortran    = ['gfortran']

charm_path = '/home/ux452912/charm/620/gnu/net/charm'

papi_path  = '/home/ux452912'
hdf5_path  = '/home/ux452912'
