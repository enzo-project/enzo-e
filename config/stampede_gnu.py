import os
from _common_search_paths import charm_path_search, grackle_path_search
is_arch_valid = 1

#python_lt_27 = 1

flags_arch = '-O3 -Wall -g'
#flags_arch = '-Wall -g'
flags_link  = '-rdynamic'

#optional fortran flag
flags_arch_fortran = '-ffixed-line-length-132'

cc   = 'gcc'
f90  = 'gfortran'

#flags_prec_single = '-fdefault-real-4 -fdefault-double-8'
flags_prec_single = ''
flags_prec_double = '-fdefault-real-8 -fdefault-double-8'


libpath_fortran = ''
libs_fortran    = ['gfortran']

home = os.environ['HOME']

hdf5_path = os.getenv('HDF5HOME',None)
if hdf5_path is not None:
	hdf5_inc = hdf5_path + '/include'
	hdf5_lib = hdf5_path + '/lib'
else:
	# the following environment variables are set by the hdf5 module
	hdf5_inc = os.environ['TACC_HDF5_INC']
	hdf5_lib = os.environ['TACC_HDF5_LIB']

boost_path = os.environ['BOOST_ROOT']
boost_inc = boost_path + '/include'
boost_lib = boost_path + '/lib'

#--------------------------------------------------
# CHARM
#
# Change charm_path below to match where your copy is.  To compile
# Charm++ on Stampede with GNU compilers, use the following:
#
#    ./build charm++ ofi-linux-x86_64 -j8  --with-production --enable-tracing
#
#--------------------------------------------------

if os.path.isdir(home + '/Charm/682/gnu/omni/charm'):
	charm_path = home + '/Charm/682/gnu/omni/charm'
else:
	charm_path = charm_path_search(home)

if ((os.getenv("TACC_PAPI_LIB", None) is not None) and
    (os.getenv("TACC_PAPI_INC", None) is not None)):
	papi_inc = os.environ["TACC_PAPI_INC"]
	papi_lin = os.environ["TACC_PAPI_LIB"]
else:
	papi_inc = home + '/include'
	papi_lib = home + '/lib'

png_path     = '/usr/lib64'

if os.path.isdir(home + '/public/Grackle/src/clib'):
	grackle_path = home + '/public/Grackle/src/clib'
else:
	grackle_path = grackle_path_search(home)
