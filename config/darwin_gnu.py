import os

is_arch_valid = 1

flags_arch = '-O3 -g'
#flags_arch = '-Wall -pg -fprofile-arcs -ftest-coverage'
#flags_arch = '-Wall -pg'
#flags_arch = '-Wall -g'
flags_link_charm = ' -rdynamic' # required for backtraces

cc  = 'gcc-mp-4.9'
f90 = 'gfortran-mp-4.9'

flags_prec_single = '-fdefault-real-4 -fdefault-double-8'
flags_prec_double = '-fdefault-real-8 -fdefault-double-8'

libpath_fortran = '.'
libs_fortran    = ['gfortran']

home = os.environ['HOME']

charm_path   = home + '/codes/charm-6.5.1'
papi_path    = '/usr/local'
hdf5_inc    = '/opt/local/include'
hdf5_lib    = '/opt/local/lib'
png_path     = '/opt/local'
grackle_path = home + '/Software/Grackle/src/clib'

