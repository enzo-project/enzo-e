import os
import sys

charm_path = '/home/bordner/charm/charm'

# CONFIGURATION

use_papi  = 0
use_hdf5  = 1
use_adios = 0

# DEFINES

define_papi  = ['CONFIG_USE_PAPI'];
define_hdf5  = ['CONFIG_USE_HDF5','H5_USE_16_API'];
define_adios = ['CONFIG_USE_ADIOS'];
define_png   = ['NO_FREETYPE'];
define_mpi   = ['CONFIG_USE_MPI'];
define_charm = ['CONFIG_USE_CHARM']
define_single = ['CONFIG_PRECISION_SINGLE']
define_double = ['CONFIG_PRECISION_DOUBLE']

# (define*_string used for IBM)

defines_string = ' -D' + define_png[0]
defines        = define_png;

if (use_papi != 0): 
	defines = defines + define_papi
	defines_string = defines_string + ' -D' + define_papi[0]

if (use_hdf5 != 0):
	defines = defines + define_hdf5
	defines_string = defines_string + ' -D' + define_hdf5[0]

if (use_adios != 0):
	defines = defines + define_adios
	defines_string = defines_string + ' -D' + define_adios[0]

#-----------------------------------------------------------------------
# PARSE ARGUMENTS
#-----------------------------------------------------------------------

arch = ARGUMENTS.get('arch','unknown')
type = ARGUMENTS.get('type','unknown')
prec = ARGUMENTS.get('prec','unknown')

if (arch == 'unknown' and "CELLO_ARCH" in os.environ):
     arch = os.environ["CELLO_ARCH"]
if (type == 'unknown' and "CELLO_TYPE" in os.environ):
     type = os.environ["CELLO_TYPE"]
if (prec == 'unknown' and "CELLO_PREC" in os.environ):
     prec = os.environ["CELLO_PREC"]

platform = arch + '-' + type

if (prec == 'unknown'): prec = 'single'

if (prec == 'single'):
	defines = defines + define_single
	defines_string = defines_string + ' -D' + define_single[0]
elif (prec == 'double'):
	defines = defines + define_double
	defines_string = defines_string + ' -D' + define_double[0]

#==================================================
# Initialize environment according to platform
#==================================================

platform_list = [
	      ['linux','mpi'],
	      ['linux','charm'],
	      ['linux','charm-perf'],
	      ['linux','serial'],
	      ['ncsa-bd','charm'],
	      ['ncsa-bd','serial'],
	      ['sdsc-triton','charm'],
	      ['sdsc-triton','serial'],
	    ]

ok = 0
for check in platform_list:
    if (platform == check[0] + '-' + check[1]):
       ok = 1

if (ok == 0):
   print "**********************************************************************"
   print
   print "Platform (<arch>-<type>) '",platform,"' is not recognized.  To specify the platform, either:"
   print
   print "1) Set the 'CELLO_ARCH' and 'CELLO_TYPE' environment variables,"
   print
   print "   or"
   print
   print "2) Use the 'arch=<arch>' and 'type=<type>' scons arguments"
   print
   print "Recognized <arch> and <type> combinations are:"
   print
   print "       ARCH           TYPE "
   print "   ------------   ------------"
   for p in platform_list:
      print "   %(arch)-12s   %(type)-12s" % \
              {'arch': p[0],     'type':p[1]}
   print
   print "**********************************************************************"
   print
   sys.exit()
       

#--------------------------------------------------
if (platform == 'linux-serial'):
#--------------------------------------------------

   flags_opt = '-g'
   flags_prec = '-m128bit-long-double'
   flags_warn = '-Wall'
   flags_cxx = flags_opt + ' ' + flags_prec + ' ' + flags_warn

   parallel_run = ""
   parallel_type = "serial"
   serial_run   = ""

   if (use_papi):
      include_path = ['#/include','/usr/local/include'];
      lib_path     = ['#/lib','/usr/local/lib'];
   else:
      include_path = ['#/include'];
      lib_path     = ['#/lib'];

   env = Environment (
      CC          = 'gcc',	
      CPPDEFINES  = defines,
      CPPPATH     = include_path,
      CXXFLAGS    = flags_cxx,
      CFLAGS      = flags_cxx,
      CXX         = 'g++',	
      ENV         = os.environ,
      FORTRAN     = 'gfortran',
      FORTRANFLAGS = flags_cxx,
      FORTRANLIBS = 'gfortran',
      FORTRANPATH = include_path,
      LIBPATH     = lib_path,
      LINKFLAGS   = flags_cxx
   )
#--------------------------------------------------
elif (platform == 'linux-mpi'):
#--------------------------------------------------
   parallel_run = "mpirun -np 4"
   serial_run   = ""
   parallel_type = "mpi"
   env = Environment (
      CC          = 'mpicc',	
      CPPDEFINES  = defines + define_mpi,
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
      CPPDEFINES  = defines + define_mpi,
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
elif (platform == 'linux-charm' or platform == 'linux-charm-perf'):
#--------------------------------------------------

   flags_opt = '-g -O3 '

#   flags_debug = '-memory charmdebug'
   flags_debug = ''

   if (platform == 'linux-charm-perf'):
	flags_charm = '-tracemode projections'
   else:
	flags_charm = ''

   flags = flags_opt + ' ' + flags_debug + ' -m128bit-long-double '
   parallel_run = charm_path + "/bin/charmrun +p4 "
   parallel_type = "charm"
   serial_run   = ""
  
   if (use_papi):
      include_path = ['#/include','/usr/local/include'];
      lib_path     = ['#/lib','/usr/local/lib'];
   else:
      include_path = ['#/include'];
      lib_path     = ['#/lib'];

   env = Environment(
      CC          = charm_path + '/bin/charmc -language charm++ '+flags+flags_charm,
      CPPDEFINES  = defines + define_charm,
      CPPPATH     = include_path,
      CXX         = charm_path + '/bin/charmc -language charm++ '+flags+flags_charm,
      CXXFLAGS    = flags,
      CFLAGS      = flags,
      ENV         = os.environ,
      FORTRAN     = 'gfortran',
      FORTRANFLAGS = flags,
      FORTRANLIBS = 'gfortran',
      FORTRANPATH = '#/include',
      LIBFLAGS     = flags,
      LIBPATH     = lib_path,
      )
   charm_builder = Builder (action="${CXX} $SOURCE; mv ${ARG}.*.h include")
   env.Append(BUILDERS = { 'CharmBuilder' : charm_builder })

#--------------------------------------------------
elif (platform == 'sdsc-triton-charm' or platform == 'sdsc-triton-serial'):
#--------------------------------------------------

   if (platform == 'sdsc-triton-charm'):
      parallel_run = charm_path + "/bin/charmrun +p4 "
      parallel_type = "charm"
   elif (platform == 'sdsc-triton-serial'):
      parallel_run = ""
      parallel_type = "serial"

   serial_run   = ""
   path_hdf5 = '/opt/hdf5/pgi'
   env = Environment (
      CC      = 'pgcc',	
      CPPDEFINES = defines + define_charm,
      CPPPATH = ['#/include', path_hdf5 + '/include'],
      CXXFLAGS = '-g -DH5_USE_16_API',
      CFLAGS   = '-g -DH5_USE_16_API',
      CXX     = 'pgCC',	
      FORTRAN = 'pgf77',
      FORTRANPATH = '#/include',
      LIBPATH = ['#/lib', path_hdf5 + '/lib'],
      LINKFLAGS = '-pgf90libs',
      ENV = { 'PATH' : os.environ['PATH'], 
              'LM_LICENSE_FILE' : os.environ['LM_LICENSE_FILE']}
   )

#--------------------------------------------------
elif (platform == 'ncsa-bd-charm' or platform == 'ncsa-bd-serial'):
#--------------------------------------------------

   if (platform == 'ncsa-bd-charm'):
      parallel_run = charm_path + "/bin/charmrun +p4 "
      parallel_type = "charm"
   elif (platform == 'ncsa-bd-serial'):
      parallel_run = ""
      parallel_type = "serial"

   serial_run    = ""

   flags_opt = '-O3 -qhot -q64'
   flags_debug = ''
   flags_fort = '-qextname -I include'
   flags = flags_opt + ' ' + flags_debug + ' '

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

   # PAPI

   path_papi = '/opt/usersoft/papi/4.1.0'

   if (use_papi != 0):
      lib_papi = path_papi + '/lib64'
      inc_papi = path_papi + '/include'
   else:
      lib_papi = ''
      inc_papi = ''
      
   
   # DEFINES 
   # (handled differently since xlC expects -Ddoh but xlf expects -D doh)


   env = Environment (
      ARFLAGS  = 'r',
      CCFLAGS = flags,
      CC      = cc + defines_string,
      CPPPATH = ['/home/bordner/include', '#/include', inc_hdf5, inc_papi],
      CXX     = cxx + defines_string,	
      CXXFLAGS = flags,
      ENV         = os.environ,
      DEFINES = '',
      FORTRANFLAGS = flags + flags_fort,
      FORTRANLIBS = ['xlf90','xlfmath','xl'],
      FORTRAN = fc,
      LIBPATH = ['#/lib','/home/bordner/lib',lib_fc,lib_hdf5,lib_papi],
      LINKFLAGS  = flags
   )

Export('env')
Export('platform')

Export('parallel_type')
Export('parallel_run')
Export('serial_run')

Export('use_papi')
Export('use_hdf5')

SConscript('src/SConscript')
SConscript('test/SConscript')

# Build tarball

# creates cello-#.#.#.tar.gz with (bin include lib)

# env = Environment(tools=['default', 'packaging'])
# title = 'Enzo-P / Cello Extreme AMR Astrophysics and Cosmology'
# env.Package( NAME           = 'cello',
#              VERSION        = '0.1.0',
#              PACKAGEVERSION = 0,
#              PACKAGETYPE    = 'targz',
#              LICENSE        = 'New BSD',
#              SUMMARY        = title,
#              DESCRIPTION    = title
#         )

