f90 = {}
cxx = {}
cc  = {}

is_arch_valid = 1

flags_arch       = '-g -pg -O3 -Wall'

# -lpthread: not needed?
# -rdynamic: required for backtraces

balancer = 'RotateLB'

flags_cxx_charm  = '-balancer ' + balancer
flags_link_charm = '-rdynamic -module ' + balancer

cc   = 'gcc'
f90  = 'gfortran'

libpath_fortran = ''
libs_fortran    = ['gfortran']

home = os.environ['HOME']

charm_path  = home + '/Charm/charm'
papi_path   = '/usr/local'
hdf5_path   = '/usr'
