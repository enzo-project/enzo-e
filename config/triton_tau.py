import os

f90 = {}
cxx = {}
cc  = {}

is_arch_valid = 1

flags_arch = '-g'
flags_link  = ''

cc   = 'icc'
f90  = 'ifort'

libpath_fortran = ''
libs_fortran    = ['imf','ifcore','ifport','stdc++']

home = os.environ['HOME']

charm_path = home + '/Charm/charm'
papi_path = ''
hdf5_path = '/opt/hdf5/intel'


#---------
# Above is triton_intel.py
# Below are updates for TAU
#---------

TAU_PATH = home + '/x86_64/bin'
cc['mpi']     = TAU_PATH + 'tau_cc.sh'
cxx['mpi']    = TAU_PATH + 'tau_cxx.sh'
