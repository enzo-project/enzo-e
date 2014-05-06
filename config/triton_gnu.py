import os 

f90 = {}
cxx = {}
cc  = {}

is_arch_valid = 1

# flags_arch = '-g -Wall'
flags_arch = '-O3'
flags_link  = '-rdynamic'

cc   = 'gcc'
f90  = 'gfortran'

libpath_fortran = ''
libs_fortran    = ['gfortran']

home = os.environ['HOME']

charm_path = home + '/Charm/charm'

papi_path  = ''
hdf5_path  = '/opt/hdf5/gnu'
