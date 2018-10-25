import os

is_arch_valid = 1

#python_lt_27 = 1

#flags_arch = '-g'
#flags_arch = '-Ofast' # ERROR: ipo linking
flags_arch = '-O3 -Wall -g'
flags_link  = '-rdynamic'

cc   = 'gcc'
f90  = 'ifort'

flags_prec_single = '-real-size 32'
flags_prec_double = '-real-size 64'

libpath_fortran = ''
libs_fortran    = ['irc', 'imf','ifcore','ifport','stdc++','intlc','svml']

home = os.environ["HOME"]
charm_path = '/home/ux452912/Charm/682/intel/net/charm'
papi_inc = home + '/include'
papi_lib = home + '/lib'
hdf5_path  = home
hdf5_inc = hdf5_path + '/include'
hdf5_lib = hdf5_path + '/lib'

png_path   = '/usr/lib64'
grackle_path = home + '/Grackle/src/clib'

