#
# Enzo-P/Cello Configuration file for NCSA Blue Waters
#
# See README.ncsa_bw for additional information

import os

f90 = {}
cxx = {}
cc  = {}

# need to load python module to get >= 2.7
# python_lt_27 = 1

node_size = 32 # for BW integer cores

is_arch_valid = 1

flags_arch = '-Wall -O3 -g -ffast-math -funroll-loops -fPIC'

flags_prec_single = ''
flags_prec_double = '-fdefault-real-8 -fdefault-double-8'

cc   = 'cc'
f90  = 'ftn'

libpath_fortran = ''
libs_fortran    = ['gfortran']

home = os.environ['HOME']

charm_path = '/sw/bw/charm/charm-v6.8.2/'
# home + '/Charm/charm.682-gni'
png_path   = '/sw/EasyBuild/software/libpng/1.6.23-CrayGNU-2016.04/'
#charm_path = os.environ['CHARM_HOME']
#png_path   = home
use_papi = 1
papi_inc="/opt/cray/papi/default/include"
papi_lib="/opt/cray/papi/default/lib64"

hdf5_path  = os.environ["HDF5_DIR"]
hdf5_inc = hdf5_path + '/include'
hdf5_lib = hdf5_path + '/lib'

grackle_path = home + '/local'

if (type == "mpi"):
   parallel_run = "aprun -n 8"
boost_path = os.environ["BOOST_ROOT"]
boost_inc = boost_path + '/include'
boost_lib = boost_path + '/lib'

node = 1
node = os.environ["CHARM_NODE"]
smp = os.environ["CHARM_SMP"]

serial_run   = 'aprun -n 1 '

if (smp == 0):
    parallel_arg = ''
    if (node == '1'):
        parallel_run = 'aprun -n 4 -N 4 -d 8 '
    elif (node == '2'):
        parallel_run = 'aprun -n 4 -N 2 -d 4 '
    elif (node == '4'):
        parallel_run = 'aprun -n 4 -N 1 -d 2 '
    else:
        parallel_run = 'error 1'
else:
    if (node == '1'):
        parallel_run = 'aprun -n 1 -N 1 -d 8 '
        parallel_arg = '++ppn 4 +setcpuaffinity '
    elif (node == '2'):
        parallel_run = 'aprun -n 2 -N 2 -d 4 '
        parallel_arg = '++ppn 2 +setcpuaffinity '
    elif (node == '4'):
        parallel_run = 'aprun -n 4 -N 4 -d 2 '
        parallel_arg = '++ppn 1 +setcpuaffinity '
    else:
        parallel_run = 'error 2'


