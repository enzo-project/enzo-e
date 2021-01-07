import os

is_arch_valid = 1

#python_lt_27 = 1

flags_arch = '-xCORE-AVX512  -O3'
#flags_arch = '-Wall -g'
flags_link  = '-rdynamic'

cc   = 'icc'
f90  = 'ifort'

flags_prec_single = ''
flags_prec_double = '-real-size 64 -double-size 64'
flags_fc = '-nofor-main'

libpath_fortran = ''
libs_fortran    = ['irc', 'imf','ifcore','ifport','stdc++','intlc','svml']

home = os.environ['HOME']

hdf5_path    = os.environ['HDF5_HOME']
hdf5_inc = hdf5_path + '/include'
hdf5_lib = hdf5_path + '/lib'

boost_path = os.environ['BOOST_HOME']
boost_inc = boost_path 
boost_lib = boost_path 

charm_path = os.environ['CHARM_HOME']

papi_inc = home + '/include'
papi_lib = home + '/lib'

png_path     = '/usr/lib64'
grackle_path = home + '/local'

node_size = 56
