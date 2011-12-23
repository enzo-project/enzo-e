     is_arch_valid = 1

     flags_arch = '-g -Wall'
     flags_link  = '-rdynamic'

     # Requires modules gnu, mpich_mx

     cc['mpi']     = 'mpicc'
     cc['serial']  = 'gcc'
     cxx['mpi']    = 'mpic++'
     cxx['serial'] = 'g++'
     f90['charm']  = 'gfortran'
     f90['mpi']    = 'gfortran'
     f90['serial'] = 'gfortran'

     libpath_fortran = ''
     libs_fortran    = ['gfortran']

     charm_path = '/home/jobordner/public/charm/charm-' + mpi_type + '-gnu'
     papi_path  = ''
     hdf5_path  = '/opt/hdf5/gnu'
