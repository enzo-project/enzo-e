import os

f90 = {}
cxx = {}
cc  = {}

is_arch_valid = 1

flags_arch = '-O3'
flags_link  = '-O3'

cc   = 'cc'
f90  = 'ftn'

libpath_fortran = ''
libs_fortran    = ['gfortran']

home_path = os.environ["HOME"]

charm_path = home + '/Charm/charm'
papi_path  = home
png_path   = home
hdf5_path  = os.environ["CRAY_HDF5_DIR"]

if (type == "mpi"):
   parallel_run = "aprun -n 8"
