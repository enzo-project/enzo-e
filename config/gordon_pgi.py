import os

f90 = {}
cxx = {}
cc  = {}

is_arch_valid = 1

flags_arch = '-O3'
flags_link = '-pgf90libs'

cc   = 'pgcc'
f90  = 'pgf90'

libpath_fortran = ''
libs_fortran    = []


charm_path = $home + '/Charm/charm'

papi_path   = home
hdf5_path   = '/opt/hdf5/pgi'
