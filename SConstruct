import os
import sys

#----------------------------------------------------------------------
# USER CONFIGURATION
#----------------------------------------------------------------------

# (TEMPORARY) Whether to bypass Block refresh and jump to output
# (for debugging adapt)

skip_refresh = 0

# Whether to print out messages with the TRACE() series of statements

trace = 0

# Whether to enable displaying messages with the DEBUG() series of statements
# Also writes messages to out.debug.<P> where P is the (physical) process rank
# Still requires the "DEBUG" group to be enabled in Monitor (that is
# Monitor::is_active("DEBUG") must be true for any output)

debug = 0

# Whether to periodically print all field values.  See
# src/Field/field_FieldBlock.cpp

debug_verbose = 0

# Whether to track dynamic memory statistics.  Can be useful, but can
# cause problems on some systems that also override new [] () / delete [] ()

memory = 1

# Enable charm++ dynamic load balancing

balance = 0

balancer = 'RotateLB'  # For testing only

# Whether to compile with -pg to use gprof for performance profiling

use_gprof = 0

# Whether to run the test programs using valgrind to check for memory leaks

use_valgrind = 0

# Whether to use Cello Performance class for collecting performance data
# (currently requires global reductions, and may not be fully function)
# (basic time data on root processor is still output)

use_performance = 0

# Whether to compile the CHARM++ version for use with the Projections
# performance tool.

use_projections = 0

# Triton MPI type (openmpi or mpich2)

mpi_type = 'mpich2'

# How many processors to run unit tests with in CHARM (MPI must be 8)

ip_charm = '4'
ip_mpi   = '8'

#----------------------------------------------------------------------
# AUTO CONFIGURATION
#----------------------------------------------------------------------

# Whether the system has the PAPI performance API installed

use_papi = 0

env = Environment()

if not env.GetOption('clean'):

     configure = Configure(env)

     if not configure.CheckCHeader('papi.h'):
          print 'PAPI not installed'
          use_papi = 0
     else:
          print 'PAPI installed'
          use_papi = 1

     env = configure.Finish()

# use_papi = 1

#-----------------------------------------------------------------------
# COMMAND-LINE ARGUMENTS
#-----------------------------------------------------------------------

# scons command line (overrides CELLO_* environment variables)

arch = ARGUMENTS.get('arch','unknown')
type = ARGUMENTS.get('type','unknown')
prec = ARGUMENTS.get('prec','unknown')

# use environment variable if scons command line not provided

if (arch == 'unknown' and "CELLO_ARCH" in os.environ):
     arch = os.environ["CELLO_ARCH"]
if (type == 'unknown' and "CELLO_TYPE" in os.environ):
     type = os.environ["CELLO_TYPE"]
if (prec == 'unknown' and "CELLO_PREC" in os.environ):
     prec = os.environ["CELLO_PREC"]

print 
print "    CELLO_ARCH scons arch=",arch
print "    CELLO_TYPE scons type=",type
print "    CELLO_PREC scons prec=",prec
print 

#----------------------------------------------------------------------
# CONFIGURATION DEFINES
#----------------------------------------------------------------------

define = {}

# Temporary defines

define_skip_refresh = ['TEMP_SKIP_REFRESH']

# Parallel type defines

define["serial"] =    []
define["mpi"]    =    ['CONFIG_USE_MPI']
define["charm"]  =    ['CONFIG_USE_CHARM']

# Precision defines

define["single"] =    ['CONFIG_PRECISION_SINGLE']
define["double"] =    ['CONFIG_PRECISION_DOUBLE']

# Performance defines

define_memory =       ['CONFIG_USE_MEMORY']
define_projections =  ['CONFIG_USE_PROJECTIONS']
define_performance =  ['CONFIG_USE_PERFORMANCE']
define_papi  =        ['CONFIG_USE_PAPI','PAPI3']

# Debugging defines

define_trace =        ['CELLO_TRACE']
define_debug =        ['CELLO_DEBUG']
define_debug_verbose = ['CELLO_DEBUG_VERBOSE']

# Library defines

define_hdf5  =        ['H5_USE_16_API']
define_png   =        ['NO_FREETYPE']

#----------------------------------------------------------------------
# ASSEMBLE DEFINES
#----------------------------------------------------------------------

defines     = []

# Parallel type configuration

if (type == 'serial' or type == 'mpi' or type == 'charm'):

     defines = defines + define[type]

else:
     print "Unrecognized parallel type ",type
     print
     print "Valid types are 'serial', 'mpi', and 'charm'"
     print
     print "The type is set using the environment variable $CELLO_TYPE"
     print "or by using 'scons type=<type>"
     sys.exit(1)

# Precision configuration

if (prec == 'single' or prec == 'double'):
     defines = defines + define[prec]
else:
     print "Unrecognized precision ",prec
     print
     print "Valid precisions are 'single' and 'double'"
     print
     print "The precision is set using the environment variable $CELLO_PREC"
     print "or by using 'scons prec=<precision>"
     sys.exit(1)

defines = defines + define_hdf5
defines = defines + define_png

charm_perf = ''

if (use_projections == 1):
     defines = defines + define_projections
     charm_perf = '-tracemode projections'

if (use_performance == 1):
     defines = defines + define_performance

flags_arch_cpp = ''
flags_config = ''
flags_cxx = ''
flags_cc = ''
flags_fc = ''
flags_link = ''
flags_cxx_charm = ''
flags_cc_charm = ''
flags_fc_charm = ''
flags_link_charm = ''

if (use_gprof == 1):
     flags_config = flags_config + ' -pg'
  
if (use_papi != 0):      defines = defines + define_papi
if (trace != 0):         defines = defines + define_trace
if (skip_refresh != 0):  defines = defines + define_skip_refresh
if (debug != 0):         defines = defines + define_debug
if (debug_verbose != 0): defines = defines + define_debug_verbose
if (memory != 0):        defines = defines + define_memory

#======================================================================
# ARCHITECTURE SETTINGS
#======================================================================

is_arch_valid = 0
sys.path.append("./config");

if   (arch == "linux-gnu"):    from linux_gnu    import *
elif (arch == "linux-gprof"):  from linux_gprof  import *
elif (arch == "linux-mpe"):    from linux_mpe    import *
elif (arch == "linux-tau"):    from linux_tau    import *
elif (arch == "ncsa-bw"):      from ncsa_bw      import *
elif (arch == "triton-gnu"):   from triton_gnu   import *
elif (arch == "triton-intel"): from triton_intel import *
elif (arch == "triton-mpe"):   from triton_mpe   import *
elif (arch == "triton-pgi"):   from triton_pgi   import *
elif (arch == "triton-tau"):   from triton_tau   import *
elif (arch == "gordon-gnu"):   from gordon_gnu   import *
elif (arch == "gordon-pgi"):   from gordon_pgi   import *
elif (arch == "gordon-intel"): from gordon_intel   import *

#======================================================================
# END ARCHITECTURE SETTINGS
#======================================================================

if (not is_arch_valid):
   print "Unrecognized architecture ",arch
   sys.exit(1)

#======================================================================
# FINAL CHARM SETUP
#======================================================================

charmc = charm_path + '/bin/charmc -language charm++ '

cxx['charm']  = charmc + charm_perf + ' '
# cc ['charm']  = charmc + charm_perf + ' '

if (balance == 1):
     flags_cxx_charm = flags_cxx_charm + ' -balancer ' + balancer
     flags_link_charm = flags_link_charm + ' -module ' + balancer

#======================================================================
# UNIT TEST SETTINGS
#======================================================================

parallel_run_args = ""

if (type == "serial"):
     serial_run   = ""
     parallel_run = ""
elif (type == "mpi"):
     serial_run   = ""
     parallel_run = "mpirun -np " + ip_mpi
elif (type == "charm"):
     serial_run   = ""
     parallel_run = charm_path + "/bin/charmrun +p" + ip_charm
     if (balance):  parallel_run_args = "+balancer " + balancer

if (use_valgrind):
     valgrind = "valgrind --leak-check=full"
     serial_run   = valgrind + " " + serial_run
     parallel_run = parallel_run + " " + valgrind

#======================================================================
# CELLO PATHS
#======================================================================

config = type

bin_path = '#/bin/'+config
lib_path = '#/lib/'+config
inc_path = '#/include'
test_path= 'test/'+config

Export('bin_path')
Export('lib_path')
Export('inc_path')
Export('test_path')
Export('ip_charm')
Export('ip_mpi')


cpppath     = [inc_path]
fortranpath = [inc_path]
libpath     = [lib_path]

#----------------------------------------------------------------------
# PAPI PATHS
#----------------------------------------------------------------------

if (use_papi):
     cpppath = cpppath + [papi_path + '/include']
     libpath = libpath + [papi_path + '/lib']

#----------------------------------------------------------------------
# HDF5 PATHS
#----------------------------------------------------------------------

cpppath = cpppath + [hdf5_path + '/include']
libpath = libpath + [hdf5_path + '/lib']

#----------------------------------------------------------------------
# FORTRAN LINK PATH
#----------------------------------------------------------------------

libpath = libpath + [libpath_fortran]

# set the Cello binary and library paths

#======================================================================
# ENVIRONMENT
#======================================================================

environ  = os.environ

cxxflags = flags_arch + ' ' + flags_arch_cpp
cxxflags = cxxflags + ' ' + flags_cxx
cxxflags = cxxflags + ' ' + flags_config
if (type=="charm"): cxxflags = cxxflags + ' ' + flags_cxx_charm

cflags   = flags_arch
cflags   = cflags + ' ' + flags_cc
cflags   = cflags + ' ' + flags_config
if (type=="charm"):cflags   = cflags + ' ' + flags_cc_charm

fortranflags = flags_arch
fortranflags = fortranflags + ' ' + flags_fc
fortranflags = fortranflags + ' ' + flags_config
if (type=="charm"):fortranflags = fortranflags + ' ' + flags_fc_charm

linkflags    = flags_arch + ' ' + flags_arch_cpp
linkflags    = linkflags + ' ' + flags_link
linkflags    = linkflags + ' ' + flags_config
if (type=="charm"):linkflags    = linkflags + ' ' + flags_link_charm

if not os.path.exists("include"):
     os.makedirs("include")
cello_def = open ("include/auto_config.def", "w")

cello_def.write ("#define CELLO_ARCH \""+arch+"\"\n")
cello_def.write ("#define CELLO_TYPE \""+type+"\"\n")
cello_def.write ("#define CELLO_PREC \""+prec+"\"\n")
cello_def.write ("#define CELLO_CC           \""+cc[type]+"\"\n")	
cello_def.write ("#define CELLO_CFLAGS       \""+cflags+"\"\n")
cello_def.write ("#define CELLO_CPPDEFINES   \""+" ".join(map(str,defines))+"\"\n")
cello_def.write ("#define CELLO_CPPPATH      \""+" ".join(map(str,cpppath))+"\"\n")
cello_def.write ("#define CELLO_CXX          \""+cxx[type]+"\"\n")	
cello_def.write ("#define CELLO_CXXFLAGS     \""+cxxflags+"\"\n")
cello_def.write ("#define CELLO_FORTRANFLAGS \""+fortranflags+"\"\n")
cello_def.write ("#define CELLO_FORTRAN      \""+f90[type]+"\"\n")
cello_def.write ("#define CELLO_FORTRANLIBS  \""+" ".join(map(str,libs_fortran))+"\"\n")
cello_def.write ("#define CELLO_FORTRANPATH  \""+" ".join(map(str,fortranpath))+"\"\n")
cello_def.write ("#define CELLO_LIBPATH      \""+" ".join(map(str,libpath))+"\"\n")
cello_def.write ("#define CELLO_LINKFLAGS    \""+linkflags+"\"\n" )

cello_def.close()

env = Environment (
     CC           = cc[type],	
     CFLAGS       = cflags,
     CPPDEFINES   = defines,
     CPPPATH      = cpppath,
     CXX          = cxx[type],	
     CXXFLAGS     = cxxflags,
     ENV          = environ,
     FORTRANFLAGS = fortranflags,
     FORTRAN      = f90[type],
     FORTRANLIBS  = libs_fortran,
     FORTRANPATH  = fortranpath,
     LIBPATH      = libpath,
     LINKFLAGS    = linkflags )

#======================================================================
# BUILDERS
#======================================================================

if (type == "charm"):
     # include files moved to include here since they are generated in
     # top-level directory
     charm_builder = Builder (action="${CXX} $SOURCE; mv ${ARG}.*.h `dirname $SOURCE`")
     env.Append(BUILDERS = { 'CharmBuilder' : charm_builder })

Export('env')
Export('type')
Export('parallel_run')
Export('parallel_run_args')
Export('serial_run')
Export('use_papi')

SConscript( 'src/SConscript',variant_dir='build/' + config )
SConscript('test/SConscript',variant_dir = test_path )

#======================================================================
# CLEANING
#======================================================================

Clean('.','test/' + type + '-' + prec)
Clean('.','bin/'  + type + '-' + prec)
Clean('.','lib/'  + type + '-' + prec)

if (type == 'charm' and use_projections == 1):
   Clean('.',Glob('bin/charm/*.projrc'))
   Clean('.',Glob('bin/charm/*.log'))
   Clean('.',Glob('bin/charm/*.sts'))
   Clean('.','charmrun')

#======================================================================
# PACKAGING
#======================================================================

# env = Environment(tools=['default', 'packaging'])
# title = 'Enzo-P / Cello Extreme AMR Astrophysics and Cosmology'
# env.Package( NAME           = 'cello',
#              VERSION        = '0.5.0',
#              PACKAGEVERSION = 0,
#              PACKAGETYPE    = 'targz',
#              LICENSE        = 'New BSD',
#              SUMMARY        = title,
#              DESCRIPTION    = title
#         )
