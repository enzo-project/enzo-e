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

cc['charm']   = 'gcc'
cc['mpi']     = 'mpicc'
cc['serial']  = 'gcc'
cxx['mpi']    = 'mpic++'
cxx['serial'] = 'g++'
f90['charm']  = 'gfortran'
f90['mpi']    = 'gfortran'
f90['serial'] = 'gfortran'

libpath_fortran = ''
libs_fortran    = ['gfortran']

charm_path  = '/home/bordner/Charm/charm'
papi_path   = '/usr/local'
hdf5_path   = '/usr'
