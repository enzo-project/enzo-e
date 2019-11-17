
import os
from _common_search_paths import charm_path_search, grackle_path_search

is_arch_valid = 1

#
#flags_arch = '-g -fprofile-arcs -ftest-coverage' # gcov
flags_arch = '-Wall -O3 -g -ffast-math -funroll-loops -fPIC'
#flags_arch = '-Wall -O3 -g'
#flags_arch = '-Wall -g'
#flags_arch = '-O3 -pg -g'
#flags_arch = '-fprofile-arcs -ftest-coverage'
#flags_arch = '-Wall -g -fsanitize=address -fno-omit-frame-pointer'
#flags_arch = '-Wall -O3 -pg'

# rdynamic required for backtraces
#flags_link_charm = '-rdynamic' 
#flags_link_charm = '-memory paranoid'
#flags_link_charm = '-fprofile-arcs' # gcov

#optional fortran flag
flags_arch_fortran = '-ffixed-line-length-132'

cc  = 'gcc '
f90 = 'gfortran'

flags_prec_single = ''
flags_prec_double = '-fdefault-real-8 -fdefault-double-8'

libpath_fortran = '/usr/lib/x86_64-linux-gnu'
libs_fortran    = ['gfortran']

home = os.getenv('HOME')

# use Charm++ with randomized message queues for debugging and stress-testing
# charm_path = home + '/Charm/charm.random'

charm_path = charm_path_search(home)

# use_papi = 1
# papi_inc = os.getenv('PAPI_INC', '/usr/include')
# papi_lib = os.getenv('PAPI_LIB', '/usr/lib')

hdf5_inc = os.getenv('HDF5_INC')
if hdf5_inc is None:
	if os.path.exists('/usr/include/hdf5.h'):
	        hdf5_inc    = '/usr/include'
	elif os.path.exists('/usr/include/hdf5/serial/hdf5.h'):
		hdf5_inc    = '/usr/include/hdf5/serial'
	else:
		raise Exception('HDF5 include file was not found.  Try setting the HDF5_INC environment variable such that $HDF5_INC/hdf5.h exists.')

hdf5_lib = os.getenv('HDF5_LIB')
if hdf5_lib is None:
	if os.path.exists('/usr/lib/libhdf5.a'):
		hdf5_lib    = '/usr/lib'
	elif os.path.exists('/usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a'):
		hdf5_lib    = '/usr/lib/x86_64-linux-gnu/hdf5/serial'
	else:
		raise Exception('HDF5 lib file was not found.  Try setting the HDF5_LIB environment variable such that $HDF5_LIB/libhdf5.a exists.')

png_path = os.getenv('LIBPNG_HOME', '/lib/x86_64-linux-gnu')

boost_inc = os.getenv('BOOST_INC', '/usr/include/boost')
boost_lib = os.getenv('BOOST_LIB', '/usr/lib/x86_64-linux-gnu')

grackle_path = grackle_path_search(home)
