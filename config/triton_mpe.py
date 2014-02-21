f90 = {}
cxx = {}
cc  = {}

is_arch_valid = 1

flags_arch = '-g'
flags_link  = ''

cc   = 'icc'
f90  = 'ifort'

libpath_fortran = '/home/jobordner/lib'
libs_fortran    = ['imf','ifcore','ifport','stdc++','lmpe', 'mpe']


charm_path = '/home/jobordner/public/charm/charm'
papi_path = ''
hdf5_path = '/opt/hdf5/intel'
