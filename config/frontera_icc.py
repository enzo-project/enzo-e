import os

is_arch_valid = 1
use_gfortran = 1
smp = 1

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

###USE GFORTRAN INSTEAD OF IFORT
if use_gfortran:
    f90 = 'gfortran'
    libpath_fortran = '/opt/apps/gcc/9.1.0/lib64'
    libs_fortran = ['gfortran']
    flags_arch_fortran = '-ffixed-line-length-132'
    flags_prec_double = '-fdefault-real-8 -fdefault-double-8'
    flags_arch = '-O3 -Wall'
    flags_fc = ''
########################

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
