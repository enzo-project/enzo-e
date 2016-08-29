import os
if "YT_DEST" in os.environ:
    is_arch_valid = 1
    home = os.environ["YT_DEST"]
else:
    is_arch_valid = 0
    home = "/no/path/"

flags_arch = '-Wall -O3 -g'
#flags_arch = '-Wall -pg -fprofile-arcs -ftest-coverage'
#flags_arch = '-Wall -pg'
#flags_arch = '-Wall -g'
flags_link_charm = ' -rdynamic' # required for backtraces

cc  = 'gcc '
f90 = 'gfortran'

flags_prec_single = '-fdefault-real-4 -fdefault-double-8'
flags_prec_double = '-fdefault-real-8 -fdefault-double-8'

libpath_fortran = '.'
libs_fortran    = ['gfortran']

charm_path   = os.path.join(home, 'src/Charm/charm')
# papi doesn't come with yt
papi_path    = '/usr/local'
hdf5_path    = home
hdf5_inc = hdf5_path + '/include'
hdf5_lib = hdf5_path + '/lib'
png_path     = home
grackle_path = home + '/src/grackle/src/clib'
