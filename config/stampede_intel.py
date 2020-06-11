import os
from _common_search_paths import charm_path_search, grackle_path_search

is_arch_valid = 1

#vec_flags = '-xCORE-AVX2 '   # mimimum instructions supported by all nodes
#vec_flags = '-xCORE-AVX512 ' # only supports SKX login & compute nodes
#vec_flags = '-xMIC-AVX512 '  # only supports KNL compute nodes
# As recommended by the stampede user-guide the following builds a "fat"
# binary including code paths for both SKX and KNL that is selected at
# execution. It's more flexible, but has (~20%) longer compile times and
# produces larger binaries and (which may be slightly slower than targetting
# one architecture.
vec_flags = '-xCORE-AVX2 -axCORE-AVX512,MIC-AVX512 '

flags_arch = vec_flags + '-Ofast -g -Wall'
flags_link = '-rdynamic'

cc   = 'icc'
f90  = 'ifort'

flags_prec_single = ''
flags_prec_double = '-real-size 64 -double-size 64'


libpath_fortran = ''
libs_fortran    = ['irc', 'imf','ifcore','ifport','stdc++','intlc','svml']

home = os.environ['HOME']

hdf5_path = os.getenv('HDF5HOME',None)
if hdf5_path is not None:
	hdf5_inc = hdf5_path + '/include'
	hdf5_lib = hdf5_path + '/lib'
else:
	# the following environment variables are set by the hdf5 module
	hdf5_inc = os.environ['TACC_HDF5_INC']
	hdf5_lib = os.environ['TACC_HDF5_LIB']

# the BOOST_ROOT environment variable is defined by the the boost module
boost_path = os.environ['BOOST_ROOT']
boost_inc = boost_path + '/include'
boost_lib = boost_path + '/lib'

#--------------------------------------------------
# CHARM
#
# Change charm_path below to match where your copy is.  To compile
# Charm++ on Stampede with Intel compilers, use the following:
#
#    ./build charm++ ofi-linux-x86_64 icc   ifort  -j8  --with-production --enable-tracing
#
#--------------------------------------------------

if os.path.isdir(home + '/Charm/682/intel/omni/charm'):
	charm_path = home + '/Charm/682/intel/omni/charm'
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
