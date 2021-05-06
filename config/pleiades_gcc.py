import os
from _common_search_paths import charm_path_search, grackle_path_search

print("##########################################################")
print("NASA's HECC Intel machiness (last verified 2021-05-05)")
print("Load the following modules:")
print("$ module load gcc/8.2 boost/1.62 mpi-hpe/mpt python3/3.5.2")
print("##########################################################")

is_arch_valid = 1

# As pointed out by the user guide, AVX512 for the Skylake and above nodes
# may degrade performance. Thus, choosing AVX2 instructions here is the
# conservative options and also ensure that the executable runs on
# Haswell, Broadwell, Skylake, and Cascadelake nodes.
flags_arch = '-O3 -g -Wall -mavx2'

#optional fortran flag
flags_arch_fortran = '-ffixed-line-length-132'

cc   = 'gcc'
f90  = 'gfortran'

flags_prec_single = ''
flags_prec_double = '-fdefault-real-8 -fdefault-double-8'


libpath_fortran = ''
libs_fortran    = ['gfortran']

home = os.environ['HOME']

hdf5_path = home + '/src/hdf5-1.10.7/install'
hdf5_inc = hdf5_path + '/include'
hdf5_lib = hdf5_path + '/lib'

# the BOOST_ROOT environment variable is defined by the the boost module
boost_path = os.environ['BOOST_ROOT']
boost_inc = boost_path + '/include'
boost_lib = boost_path + '/lib'

#--------------------------------------------------
# CHARM
#
# Change charm_path below to match where your copy is.  To compile
# Charm++ on Pleiades with Intel compilers, use the following:
#
#    ./build charm++ mpi-linux-x86_64 smp gcc gfortran --with-production -j8
#
#--------------------------------------------------

if os.path.isdir(home + '/src/charm-v6.10.2'):
	charm_path = home + '/src/charm-v6.10.2'
else:
	charm_path = charm_path_search(home)

smp = 1

## no PAPI for now
# papi_inc = home + '/include'
# papi_lib = home + '/lib'

png_path     = '/usr/lib64'

## to enable Grackle set following value to 1 and set path below
use_grackle = 0
if os.path.isdir(home + '/public/Grackle/src/clib'):
	grackle_path = home + '/public/Grackle/src/clib'
else:
	grackle_path = grackle_path_search(home)
