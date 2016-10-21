
import os

f90 = {}
cxx = {}
cc  = {}

python_lt_27 = 1

node_size = 32 # for BW integer cores

is_arch_valid = 1

flags_arch = '-O3'
flags_link  = ''

flags_prec_single = ''
flags_prec_double = ''

cc   = 'cc'
f90  = 'ftn'

libpath_fortran = ''
libs_fortran    = ['gfortran']

home = os.environ['HOME']

charm_path = home + '/Charm/charm'
png_path   = home
#papi_path  = home
hdf5_path  = os.environ["HDF5_DIR"]
hdf5_inc = hdf5_path + '/include'
hdf5_lib = hdf5_path + '/lib'

if (type == "mpi"):
   parallel_run = "aprun -n 8"
