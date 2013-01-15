f90 = {}
cxx = {}
cc  = {}

is_arch_valid = 1

flags_arch = ''
flags_link  = ''

cc['charm']   = 'cc -h gnu'
cc['mpi']     = 'cc -h gnu'
cc['serial']  = 'cc -h gnu'
cxx['mpi']    = 'CC'
cxx['serial'] = 'CC'
f90['charm']  = 'ftn -rm'
f90['mpi']    = 'ftn -rm'
f90['serial'] = 'ftn -rm'

libpath_fortran = ''
libs_fortran    = []

charm_path = '/u/sciteam/bordner/Charm/charm'
papi_path  = '/u/sciteam/bordner'
hdf5_path  = '/opt/cray/hdf5/default/cray/74'

if (type == "mpi"):
   parallel_run = "aprun -n 8"
