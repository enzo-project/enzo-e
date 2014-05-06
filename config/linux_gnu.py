import os

f90 = {}
cxx = ''
cc  = {}

is_arch_valid = 1

flags_arch       = '-O3'
#flags_arch       = '-g -Wall'

#flags_arch_cpp = '-memory paranoid'
# -lpthread: not needed?
# -rdynamic: required for backtraces

flags_link_charm = ' -rdynamic'

cc  = 'gcc '
f90 = 'gfortran'

libpath_fortran = '.'
libs_fortran    = ['gfortran']

home = os.environ['HOME']

charm_path  = home + '/Charm/charm'
papi_path   = '/usr/local'
hdf5_path   = '/usr'
png_path    = '/lib/x86_64-linux-gnu'

