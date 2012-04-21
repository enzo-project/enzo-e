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

charm_path = '/home/jobordner/public/charm/charm-620-intel'
papi_path = ''
hdf5_path = '/opt/hdf5/intel'


#---------
# Above is triton_intel.py
# Below are updates for TAU
#---------

TAU_PATH = '/home/jobordner/x86_64/bin/'
cc['mpi']     = TAU_PATH + 'tau_cc.sh'
cxx['mpi']    = TAU_PATH + 'tau_cxx.sh'
