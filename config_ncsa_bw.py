f90 = {}
cxx = {}
cc  = {}

is_arch_valid = 1

flags_arch = ''
flags_link  = ''

cc['mpi']     = 'cc -h gnu'
cc['serial']  = 'cc -h gnu'
cxx['mpi']    = 'CC'
cxx['serial'] = 'CC'
f90['charm']  = 'ftn -rm'
f90['mpi']    = 'ftn -rm'
f90['serial'] = 'ftn -rm'

libpath_fortran = '/home/cpv/tra61/lib'
libs_fortran    = []

charm_path = '/home/cpv/tra61/charm'
papi_path  = ''
hdf5_path  = '/opt/cray/hdf5/1.8.6/cray/73'

if (type == "mpi"):
   parallel_run = "aprun -n 8"
