import os
import sys

#----------------------------------------------------------------------
# CONFIGURATION
#----------------------------------------------------------------------

use_hdf5     = 1
use_png      = 1
use_papi     = 0
use_valgrind = 0

#-----------------------------------------------------------------------
# PARSE ARGUMENTS
#-----------------------------------------------------------------------

#--------------------------------------------------
#    arch=<arch> (CELLO_ARCH)
#
#      linux         generic linux box
#      ncsa-bd       NCSA's IBM blue-drop
#      sdsc-triton   SDSC's Triton
#--------------------------------------------------
#    type=<type> (CELLO_TYPE)
#
#      serial        Compile without parallelism
#      charm         Compile with CHARM++ parallelism
#      mpi           Compile with MPI parallelism
#--------------------------------------------------
#    prec=<prec> (CELLO_PREC)
#
#      single        Compile Enzo-P for single precision
#      double        Compile Enzo-P for double precision
#--------------------------------------------------

# scons command line (overrides CELLO_* environment variables)

arch = ARGUMENTS.get('arch','unknown')
type = ARGUMENTS.get('type','unknown')
prec = ARGUMENTS.get('prec','unknown')

parallel_type = type

# environment variable (default if scons command line not provided)

if (arch == 'unknown' and "CELLO_ARCH" in os.environ):
     arch = os.environ["CELLO_ARCH"]
if (type == 'unknown' and "CELLO_TYPE" in os.environ):
     type = os.environ["CELLO_TYPE"]
if (prec == 'unknown' and "CELLO_PREC" in os.environ):
     prec = os.environ["CELLO_PREC"]

#----------------------------------------------------------------------
# DEFINES
#----------------------------------------------------------------------

define_serial =  [];
define_mpi    =  ['CONFIG_USE_MPI'];
define_charm  =  ['CONFIG_USE_CHARM']

define_single = ['CONFIG_PRECISION_SINGLE']
define_double = ['CONFIG_PRECISION_DOUBLE']

define_hdf5  =  ['CONFIG_USE_HDF5','H5_USE_16_API'];
define_png   =  ['NO_FREETYPE'];
define_papi  =  ['CONFIG_USE_PAPI'];

defines     = []
defines_xlc = ""

#--------------------------------------------------

if (type == 'serial'):
	defines     = defines
	defines_xlc = defines_xlc
elif (type == 'mpi'):
	defines     = defines             + define_mpi
	defines_xlc = defines_xlc + ' -D' + define_mpi[0]
elif (type == 'charm'):
	defines     = defines             + define_charm
	defines_xlc = defines_xlc + ' -D' + define_charm[0]
else:
	print "Unrecognized parallel type ",type
	sys.exit(1)

#--------------------------------------------------

if (prec == 'single'):
	defines = defines + define_single
	defines_xlc = defines_xlc + ' -D' + define_single[0]
elif (prec == 'double'):
	defines = defines + define_double
	defines_xlc = defines_xlc + ' -D' + define_double[0]
else:
	print "Unrecognized precision ",prec
	sys.exit(1)

#--------------------------------------------------

if (use_papi != 0): 
	defines     = defines             + define_papi
	defines_xlc = defines_xlc + ' -D' + define_papi[0]

#--------------------------------------------------

if (use_hdf5 != 0):
	defines     = defines     +         define_hdf5
	defines_xlc = defines_xlc + ' -D' + define_hdf5[0]

#--------------------------------------------------

if (use_png != 0):
	defines     = defines             + define_png;
	defines_xlc = defines_xlc + ' -D' + define_png[0]

#--------------------------------------------------


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

charm_path = '/home/bordner/charm/charm'  # arch

cpppath     = ['#/include'];
fortranpath = ['#/include'];
libpath     = ['#/lib'];

#======================================================================
# ARCHITECTURE SETTINGS
#======================================================================

if (arch == "linux"):

   fortran_serial = 'gfortran'
   fortran_mpi    = 'gfortran'
   fortran_charm  = 'gfortran'

   cxx_serial = 'g++'
   cxx_mpi    = 'mpic++'
   cxx_charm  = charm_path + '/bin/charmc -language charm++ '

   cc_serial  = 'gcc'
   cc_mpi     = 'mpicc'
   cc_charm   = charm_path + '/bin/charmc -language charm++ '

   cppdefines = defines

   papi_home = '/usr/local'

#   flags_opt = '-g -O3 -pg'
   flags_opt = '-g -O3 '
#   flags_opt  = '-g'
   flags_prec = '-m128bit-long-double'
   flags_warn = '-Wall'

# compiler-dependent flags
   cxxflags_opt  = flags_opt
   cxxflags_prec = flags_prec
   cxxflags_warn = flags_warn

   cflags_opt  = flags_opt
   cflags_prec = flags_prec
   cflags_warn = flags_warn

   fortranflags_opt  = flags_opt
   fortranflags_prec = flags_prec
   fortranflags_warn = flags_warn

   linkflags_opt  = flags_opt
   linkflags_prec = flags_prec
   linkflags_warn = flags_warn

elif (arch == "ncsa-bd"):

   print "unfinished"

elif (arch == "sdsc-triton"):

   print "unfinished"


#======================================================================
# PARALLELISM SETTINGS
#======================================================================
flags_type = ''
if (type == "serial"):
   cxx          = cxx_serial
   cc           = cc_serial
   fortran      = fortran_serial
   serial_run   = ""
   parallel_run = ""
elif (type == "mpi"):
   cxx          = cxx_mpi
   cc           = cc_mpi
   fortran      = fortran_mpi
   serial_run   = ""
   parallel_run = "mpirun -np 4"
elif (type == "charm"):
   cxx          = cxx_charm
   cc           = cc_charm
   fortran      = fortran_charm
   serial_run   = ""
   parallel_run = charm_path + "/bin/charmrun +p4 "
#   flags_type =  '-tracemode projections'

#======================================================================
# PERFORMANCE SETTINGS
#======================================================================

if (use_papi):
   cpppath     = cpppath     + (papi_home + '/include')

#======================================================================
# DEBUG SETTINGS
#======================================================================

if (use_valgrind):
   parallel_run = parallel_run + " valgrind"
   serial_run   = "valgrind " + serial_run


environ  = os.environ

cxxflags = cxxflags_opt + ' ' + cxxflags_prec + ' ' + cxxflags_warn
cflags   = cflags_opt + ' '  + cflags_prec + ' ' + cflags_warn
fortranflags = fortranflags_opt + ' ' + fortranflags_prec + ' ' + fortranflags_warn
linkflags = linkflags_opt + ' ' + linkflags_prec + ' ' + linkflags_warn

platform = arch + '-' + type

#======================================================================
# ENVIRONMENT
#======================================================================

env = Environment (
      CC          = cc,	
      CFLAGS      = cflags,
      CPPDEFINES  = cppdefines,
      CPPPATH     = cpppath,
      CXX         = cxx,	
      CXXFLAGS    = cxxflags,
      ENV         = environ,
      FORTRANFLAGS = fortranflags,
      FORTRAN     = fortran,
      FORTRANLIBS = fortran,
      FORTRANPATH = fortranpath,
      LIBPATH     = libpath,
      LINKFLAGS   = linkflags )

#======================================================================
# BUILDERS
#======================================================================

if (type == "charm"):
   charm_builder = Builder (action="${CXX} $SOURCE; mv ${ARG}.*.h include")
   env.Append(BUILDERS = { 'CharmBuilder' : charm_builder })

#--------------------------------------------------
if (platform == 'sdsc-triton-charm' or platform == 'sdsc-triton-serial'):
#--------------------------------------------------

   path_hdf5 = '/opt/hdf5/pgi'
   env = Environment (
      CC      = 'pgcc',	
      CPPDEFINES = defines,
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
   elif (platform == 'ncsa-bd-serial'):
      parallel_run = ""

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
      CC      = cc + defines_xlc,
      CPPPATH = ['/home/bordner/include', '#/include', inc_hdf5, inc_papi],
      CXX     = cxx + defines_xlc,	
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

