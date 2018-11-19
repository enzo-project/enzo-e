import os

is_arch_valid = 1

#python_lt_27 = 1

flags_arch = '-O3 -Wall -g'
#flags_arch = '-Wall -g'
flags_link  = '-rdynamic'

cc   = 'gcc'
f90  = 'gfortran'

#flags_prec_single = '-fdefault-real-4 -fdefault-double-8'
flags_prec_single = ''
flags_prec_double = '-fdefault-real-8 -fdefault-double-8'


libpath_fortran = ''
libs_fortran    = ['gfortran']

home = os.environ['HOME']

hdf5_path    = os.environ['HDF5HOME']
hdf5_inc = hdf5_path + '/include'
hdf5_lib = hdf5_path + '/lib'

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

charm_path   = home + '/Charm/682/gnu/omni/charm'
papi_inc = home + '/include'
papi_lib = home + '/lib'

png_path     = '/usr/lib64'
grackle_path = home + '/public/Grackle/src/clib'
