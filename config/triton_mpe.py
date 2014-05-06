import os

f90 = {}
cxx = {}
cc  = {}

is_arch_valid = 1

flags_arch = '-g'
flags_link  = ''

cc   = 'icc'
f90  = 'ifort'

home = os.environ['HOME']

libpath_fortran = home + '/lib'
libs_fortran    = ['imf','ifcore','ifport','stdc++','lmpe', 'mpe']


charm_path = home + '/Charm/charm'
papi_path = ''
hdf5_path = '/opt/hdf5/intel'
