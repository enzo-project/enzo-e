f90 = {}
cxx = {}
cc  = {}

is_arch_valid = 1

flags_arch = '-g'
#flags_arch = '-Ofast' # ERROR: ipo linking
#flags_arch = '-O3'
flags_link  = '-g'

cc['charm']   = 'icc'
cc['mpi']     = 'mpicc'
cc['serial']  = 'icc'
cxx['mpi']    = 'mpicxx'
cxx['serial'] = 'icpc'
f90['charm']  = 'ifort'
f90['mpi']    = 'ifort'
f90['serial'] = 'ifort'

libpath_fortran = ''
libs_fortran    = ['imf','ifcore','ifport','stdc++','intlc']

charm_path = '/home/ux452912/Charm/charm'

papi_path = '/home/ux452912'
hdf5_path = '/opt/hdf5/intel'
