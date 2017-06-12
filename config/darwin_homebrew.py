import os

is_arch_valid = 1

flags_arch = '-O3 -g'
flags_link_charm = ' -rdynamic' # required for backtraces

cc  = 'cc'
# you will need to do "brew install gfortran"
f90 = 'gfortran'

flags_prec_single = '-fdefault-real-4 -fdefault-double-8'
flags_prec_double = '-fdefault-real-8 -fdefault-double-8'

libpath_fortran = '/usr/local/lib/gcc/6'
libs_fortran    = ['gfortran']

home = os.environ['HOME']

# edit to point to path of local charm++ installation
charm_path   = home + '/Documents/charm-6.7.1'
papi_inc    = '/usr/local/include'
papi_lib    = '/usr/local/lib'
hdf5_inc    = '/usr/local/include'
hdf5_lib    = '/usr/local/lib'
png_path     = '/usr/local'

# edit to point to path of local grackle installation
grackle_path = home + '/Documents/grackle/src/clib'

