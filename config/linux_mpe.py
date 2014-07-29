import os

f90 = {}
cxx = {}
cc  = {}

is_arch_valid = 1

#flags_arch = '-g -Wall'
flags_arch = '-O3 -Wall'
flags_link = '-rdynamic -lpthread'

flags_prec_single = '-fdefault-real-4 -fdefault-double-8'
flags_prec_double = '-fdefault-real-8 -fdefault-double-8'

cc   = 'gcc'
f90  = 'gfortran'

libpath_fortran = ''
libs_fortran    = ['gfortran','lmpe','mpe']

home = os.environ['HOME']

charm_path  = home + '/Charm/charm'
papi_path   = '/usr/local'
hdf5_path   = '/usr'
