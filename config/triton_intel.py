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

charm_path = home + '/Charm/charm'

papi_path = ''
hdf5_path = '/opt/hdf5/intel'
