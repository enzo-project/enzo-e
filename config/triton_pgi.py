f90 = {}
cxx = {}
cc  = {}

is_arch_valid = 1

flags_arch = '-Ktrap=fp -g -O3'
flags_link = '-pgf90libs'

cc['charm']   = 'pgcc'
cc['mpi']     = 'mpicc'
cc['serial']  = 'pgcc'
cxx['mpi']    = 'mpicxx'
cxx['serial'] = 'pgCC'
f90['charm']  = 'pgf90'
f90['mpi']    = 'pgf90'
f90['serial'] = 'pgf90'

libpath_fortran = ''
libs_fortran    = []


charm_path  = '/home/jobordner/public/charm/charm'

papi_path   = ''
hdf5_path   = '/opt/hdf5/pgi'
