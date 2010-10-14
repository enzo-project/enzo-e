import os
import sys

TOPDIR  = GetLaunchDir()

charm_path = '/home/bordner/charm/charm'

# Initialize platform

platform='unknown'

platform = ARGUMENTS.get('platform','unknown')

if (platform == 'unknown' and "CELLO_PLATFORM" in os.environ):
     platform = os.environ["CELLO_PLATFORM"]

#==================================================
# Initialize environment according to platform
#==================================================

platform_list = [
	      'linux-ampi',
	      'linux-mpi',
	      'linux-charm',
	      'linux-charm-perf',
	      'linux-serial',
	      'ncsa-bd',
	      'ncsa-bd-serial'
	      'sdsc-triton',
	      'sdsc-triton-serial',
	    ]

ok = 0
for check in platform_list:
    if (platform == check):
       ok = 1
if (ok == 0):
   print "**********************************************************************"
   print
   print "Platform '",platform,"' is not recognized.  To specify the platform, either:"
   print
   print "1) Set the 'CELLO_PLATFORM' environment variable,"
   print
   print "   or"
   print
   print "2) Use the 'platform=<platform>' scons argument"
   print
   print "Recognized platforms are:"
   print
   for p in platform_list:
      print "   ",p
   print
   print "**********************************************************************"
   print
   sys.exit()
       

#--------------------------------------------------
if (platform == 'linux-serial'):
#--------------------------------------------------
   parallel_run = ""
   parallel_type = "serial"
   serial_run   = ""
   env = Environment (
      CC          = 'gcc',	
      CPPDEFINES = ['NO_FREETYPE'],
      CPPPATH     = '#/include',
      CXXFLAGS    = '-Wall -g  -m128bit-long-double',
      CFLAGS      = '-Wall -g  -m128bit-long-double',
      CXX         = 'g++',	
      ENV         = os.environ,
      FORTRAN     = 'gfortran',
      FORTRANLIBS = 'gfortran',
      FORTRANPATH = '#/include',
      LIBPATH     = '#/lib',
   )
#--------------------------------------------------
elif (platform == 'linux-mpi'):
#--------------------------------------------------
   parallel_run = "mpirun -np 4"
   serial_run   = ""
   parallel_type = "mpi"
   env = Environment (
      CC          = 'mpicc',	
      CPPDEFINES = ['NO_FREETYPE','CONFIG_USE_MPI'],
      CPPPATH     = '#/include',
      CXXFLAGS    = '-Wall -g  -m128bit-long-double',
      CFLAGS      = '-Wall -g  -m128bit-long-double',
      CXX         = 'mpic++',	
      ENV         = os.environ,
      FORTRAN     = 'gfortran',
      FORTRANLIBS = 'gfortran',
      FORTRANPATH = '#/include',
      LIBPATH     = '#/lib',
   )
#--------------------------------------------------
elif (platform == 'linux-mpi-valgrind'):
#--------------------------------------------------
   parallel_run = "mpirun -np 4 valgrind"
   serial_run   = "valgrind "
   parallel_type = "mpi"
   env = Environment (
      CC          = 'mpicc',	
      CPPDEFINES = ['NO_FREETYPE','CONFIG_USE_MPI'],
      CPPPATH     = '#/include',
      CXXFLAGS    = '-Wall -g  -m128bit-long-double',
      CFLAGS      = '-Wall -g  -m128bit-long-double',
      CXX         = 'mpic++',	
      FORTRAN     = 'gfortran',
      FORTRANLIBS = 'gfortran',
      FORTRANPATH = '#/include',
      LIBPATH     = '#/lib',
   )
#--------------------------------------------------
elif (platform == 'linux-ampi'):
#--------------------------------------------------
   parallel_run = charm_path + "/bin/charmrun +p4 "
   serial_run   = ""
   parallel_type = "ampi"
  
   env = Environment(
      CC          = charm_path + '/bin/charmc -language ampi',
      CPPDEFINES = ['NO_FREETYPE','CONFIG_USE_MPI'],
      CPPPATH     = '#/include',
      CXX         = charm_path + '/bin/charmc -language ampi',
      CXXFLAGS    = '-g',
      CFLAGS      = '-g',
      ENV         = os.environ,
      FORTRANFLAGS  = '-g',
      FORTRAN     = 'gfortran',
      FORTRANLIBS = 'gfortran',
      FORTRANPATH = '#/include',
      LINKFLAGS     = '-g',
      LIBPATH     = '#/lib',
   )
#--------------------------------------------------
elif (platform == 'linux-charm' or platform == 'linux-charm-perf'):
#--------------------------------------------------

   flags_opt = '-g '

#   flags_debug = '-memory charmdebug'
   flags_debug = ''

   if (platform == 'linux-charm-perf'):
	flags_charm = '-tracemode projections'
   else:
	flags_charm = ''

   flags = flags_opt + ' ' + flags_debug + ' '
   parallel_run = charm_path + "/bin/charmrun +p4 "
   parallel_type = "charm"
   serial_run   = ""
  
   env = Environment(
      CC          = charm_path + '/bin/charmc -language charm++ '+flags+flags_charm,
      CPPDEFINES = ['NO_FREETYPE','CONFIG_USE_CHARM'],
      CPPPATH     = '#/include',
      CXX         = charm_path + '/bin/charmc -language charm++ '+flags+flags_charm,
      CXXFLAGS    = flags,
      CFLAGS      = flags,
      ENV         = os.environ,
      FORTRAN     = 'gfortran',
      FORTRANFLAGS = flags,
      FORTRANLIBS = 'gfortran',
      FORTRANPATH = '#/include',
      LIBFLAGS     = flags,
      LIBPATH     = '#/lib' )
   charm_builder = Builder (action="${CXX} $SOURCE; mv ${ARG}.*.h include")
   env.Append(BUILDERS = { 'CharmBuilder' : charm_builder })
#--------------------------------------------------
elif (platform == 'sdsc-triton' or platform == 'sdsc-triton-serial'):
#--------------------------------------------------
   if (platform == 'sdsc-triton'):
      parallel_run = charm_path + "/bin/charmrun +p4 "
      parallel_type = "charm"
   elif (platform == 'sdsc-triton-serial'):
      parallel_run = ""
      parallel_type = "serial"

   serial_run   = ""
   path_hdf5 = '/opt/hdf5/pgi'
   env = Environment (
      CC      = 'pgcc',	
      CPPDEFINES = ['NO_FREETYPE'],
      CPPPATH = ['#/include', path_hdf5 + '/include'],
      CXXFLAGS = '-g -DH5_USE_16_API',
      CFLAGS   = '-g -DH5_USE_16_API',
      CXX     = 'pgCC',	
      FORTRAN = 'pgf77',
      FORTRANPATH = '#/include',
      LIBPATH = ['#/lib', path_hdf5 + '/lib'],
      LINKFLAGS = '-pgf90libs',
      ENV = {'PATH' : os.environ['PATH'], 'LM_LICENSE_FILE' : os.environ['LM_LICENSE_FILE']},
   )
#--------------------------------------------------
elif (platform == 'ncsa-bd' or platform == 'ncsa-bd-serial'):
#--------------------------------------------------

   if (platform == 'ncsa-bd'):
      parallel_run = charm_path + "/bin/charmrun +p4 "
      parallel_type = "charm"
   elif (platform == 'ncsa-bd-serial'):
      parallel_run = ""
      parallel_type = "serial"

   serial_run    = ""

#     flags_opt = '-O3 -qhot -q64'
   flags_opt = '-q64'
   flags_debug = '-g'
   flags_fort = '-qextname -I include'
   flags = flags_opt + ' ' + flags_debug + ' '

   defines = '-D NO_FREETYPE -DH5_USE_16_API '

   # Compilers

   path_fc   = '/opt/ibmcmp/xlf/13.1'
   path_cc   = '/opt/ibmcmp/vac/11.1'
   path_cxx  = '/opt/ibmcmp/vacpp/11.1'

   cc = path_cc + '/bin/xlc_r'
   fc = path_fc + '/bin/xlf_r'
   cxx = path_cxx + '/bin/xlC_r'

   lib_fc = path_fc + '/lib64'

   # HDF5 

   path_hdf5 = '/opt/hdf5-1.8.4-patch1-64bit'
   lib_hdf5 = path_hdf5 + '/lib'
   inc_hdf5 = path_hdf5 + '/include'

   env = Environment (
      ARFLAGS  = 'r',
      CCFLAGS = flags + defines,
      CC      = cc,
      CPPPATH = ['/home/bordner/include', '#/include', inc_hdf5],
      CXX     = cxx,	
      ENV         = os.environ,
      CXXFLAGS = flags + defines,
      DEFINES = '',
      FORTRANFLAGS = flags + flags_fort,
      FORTRANLIBS = ['xlf90','xlfmath'],
      FORTRAN = fc,
      LIBPATH = ['#/lib','/home/bordner/lib',lib_fc,lib_hdf5],
      LINKFLAGS  = flags
   )

# Generate the CELLO_PLATFORM file (should be a target)

f = open('CELLO_PLATFORM','w')
f.write(platform + '\n')
f.close()


Export('env')
Export('platform')
Export('parallel_type')
Export('parallel_run')
Export('serial_run')

SConscript('src/SConscript')
SConscript('test/SConscript')

# Build tarball

# creates cello-#.#.#.tar.gz with (bin include lib)

# env = Environment(tools=['default', 'packaging'])
# env.Package( NAME           = 'cello',
#              VERSION        = '0.1.0',
#              PACKAGEVERSION = 0,
#              PACKAGETYPE    = 'targz',
#              LICENSE        = 'New BSD',
#              SUMMARY        = 'Cello Extreme AMR Framework',
#              DESCRIPTION    = 'Cello Extreme AMR Framework'
#         )

