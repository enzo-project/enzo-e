f90 = {}
cxx = {}
cc  = {}

is_arch_valid = 1

# flags_arch = '-g -Wall'
flags_arch = '-O3'
flags_link  = '-rdynamic'

cc   = 'gcc'
f90  = 'gfortran'

libpath_fortran = ''
libs_fortran    = ['gfortran']

charm_path = '/home/ux452912/Charm/charm'

papi_path  = '/home/ux452912'
hdf5_path  = '/home/ux452912'
png_path   = '/usr/lib64'
