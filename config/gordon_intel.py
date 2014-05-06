import os

f90 = {}
cxx = {}
cc  = {}

is_arch_valid = 1

flags_arch = '-g'
#flags_arch = '-Ofast' # ERROR: ipo linking
#flags_arch = '-O3'
flags_link  = '-g'

cc   = 'icc'
f90  = 'ifort'

libpath_fortran = ''
libs_fortran    = ['imf','ifcore','ifport','stdc++','intlc']

charm_path = home + '/Charm/charm'

papi_path = home
hdf5_path = '/opt/hdf5/intel'
