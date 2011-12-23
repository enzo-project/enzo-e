     is_arch_valid = 1

     flags_arch = '-Ktrap=fp -g'
     flags_link = '-pgf90libs'

     # Requires modules pgi, mpich_mx

     cc['mpi']     = 'mpicc'
     cc['serial']  = 'pgcc'
     cxx['mpi']    = 'mpicxx'
     cxx['serial'] = 'pgCC'
     f90['charm']  = 'pgf90'
     f90['mpi']    = 'pgf90'
     f90['serial'] = 'pgf90'

     libpath_fortran = ''
     libs_fortran    = []

     charm_path  = '/home/jobordner/public/charm/charm-' + mpi_type + '-pgi'
     papi_path   = ''
     hdf5_path   = '/opt/hdf5/pgi'
