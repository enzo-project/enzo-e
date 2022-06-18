import os
from _common_search_paths import charm_path_search, grackle_path_search

is_arch_valid = 1
use_gfortran = 0
smp = 0

flags_arch = '-Wall -O3 -g'
#flags_arch = '-fprofile-arcs -ftest-coverage'
#flags_arch = '-Wall -g'
#flags_arch = '-Wall -g -fsanitize=address -fno-omit-frame-pointer'
#flags_arch = '-Wall -O3 -pg'

# rdynamic required for backtraces
#flags_link_charm = '-rdynamic' 
#flags_link_charm = '-memory paranoid' 

intel_dir = '/cm/shared/apps/spack/cpu/opt/spack/linux-centos8-zen/gcc-8.3.1/intel-19.1.1.217-4d42ptjd6wsnh5bgbzcv6lp44vxpjwut/compilers_and_libraries_2020.1.217/linux/bin/intel64'
cc  = intel_dir + '/icc'
f90 = intel_dir + '/ifort'

flags_prec_single = ''
flags_prec_double = '-r8'

libpath_fortran = '/cm/shared/apps/spack/cpu/opt/spack/linux-centos8-zen/gcc-8.3.1/intel-19.1.1.217-4d42ptjd6wsnh5bgbzcv6lp44vxpjwut/compilers_and_libraries_2020.1.217/linux/compiler/lib/intel64'
libs_fortran    = ['ifcore', 'ifport']


#USE GFORTRAN INSTEAD OF IFORT
if use_gfortran:
    f90 = 'gfortran'
    libpath_fortran = '/cm/shared/apps/spack/cpu/opt/spack/linux-centos8-zen/gcc-8.3.1/gcc-10.2.0-n7su7jf54rc7l2ozegds5xksy6qhrjin/lib64'
    libs_fortran = ['gfortran']
    flags_arch_fortran = '-ffixed-line-length-132'
    flags_prec_double = '-fdefault-real-8 -fdefault-double-8'
    flags_arch = '-O3 -Wall'
    flags_fc = ''

#############################

home = os.getenv('HOME')

charm_path = charm_path_search(home)

use_papi=0
papi_inc = '/usr/local/include'
papi_lib = '/usr/local/lib'

boost_path = os.getenv('BOOST_HOME')
boost_inc  = boost_path + '/include'
boost_lib  = boost_path + '/lib'

hdf5_inc = os.getenv('HDF5HOME') + '/include'
if hdf5_inc is None:
	if os.path.exists('/usr/include/hdf5.h'):
	        hdf5_inc    = '/usr/include'
	elif os.path.exists('/usr/include/hdf5/serial/hdf5.h'):
		hdf5_inc    = '/usr/include/hdf5/serial'
	else:
		raise Exception('HDF5 include file was not found.  Try setting the HDF5_INC environment variable such that $HDF5_INC/hdf5.h exists.')

hdf5_lib = os.getenv('HDF5HOME') + '/lib'
if hdf5_lib is None:
	if os.path.exists('/usr/lib/libhdf5.a'):
		hdf5_lib    = '/usr/lib'
	elif os.path.exists('/usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a'):
		hdf5_lib    = '/usr/lib/x86_64-linux-gnu/hdf5/serial'
	else:
		raise Exception('HDF5 lib file was not found.  Try setting the HDF5_LIB environment variable such that $HDF5_LIB/libhdf5.a exists.')

png_path = os.getenv('LIBPNG_HOME')
if png_path is None:
	png_path     = '/lib/x86_64-linux-gnu'

grackle_path = os.getenv('GRACKLE_HOME') 
