f90 = {}
cxx = {}
cc  = {}

is_arch_valid = 1

flags_arch = '-g'
flags_link  = ''

cc['mpi']     = 'mpicc'
cc['serial']  = 'icc'
cxx['mpi']    = 'mpicxx'
cxx['serial'] = 'icpc'
f90['charm']  = 'ifort'
f90['mpi']    = 'ifort'
f90['serial'] = 'ifort'

libpath_fortran = ''
libs_fortran    = ['imf','ifcore','ifport','stdc++']

charm_path = '/home/ux452912/charm/charm'

papi_path = '/home/ux452912'
hdf5_path = '/opt/hdf5/intel'
