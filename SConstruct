import os
import sys

#----------------------------------------------------------------------
# USER CONFIGURATION
#----------------------------------------------------------------------

# whether to compile "new initial" code (TEMPORARY)

new_initial = 1

# Whether to print out messages with the TRACE() series of statements

trace           = 0

# Whether to enable displaying messages with the DEBUG() series of statements
# Also writes messages to out.debug.<P> where P is the (physical) process rank
# Still requires the "DEBUG" group to be enabled in Monitor (that is
# Monitor::is_active("DEBUG") must be true for any output)

debug           = 1

# Whether to periodically print all field values.  See
# src/Field/field_FieldBlock.cpp

debug_verbose   = 0

# Whether to track dynamic memory statistics.  Can be useful, but can
# cause problems on some systems that also override new [] () / delete [] ()

memory          = 1

# Limit CHARM++ load balancing to only when AtSync() is called.  See
# src/Mesh/mesh_Block.cpp

atsync          = 0

# Whether to compile with -pg to use gprof for performance profiling

use_gprof       = 0

# Whether to run the test programs using valgrind to check for memory leaks

use_valgrind    = 0

# Whether to compile the CHARM++ version for use with the Projections
# performance tool.

use_projections = 0

# Triton MPI type (openmpi or mpich2)

mpi_type = 'mpich2'

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

use_papi = 0

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

# Parallel type defines

define["serial"] =        []
define["mpi"]    =        ['CONFIG_USE_MPI']
define["charm"]  =        ['CONFIG_USE_CHARM']

# Precision defines

define["single"] =        ['CONFIG_PRECISION_SINGLE']
define["double"] =        ['CONFIG_PRECISION_DOUBLE']

# Performance defines

define_atsync =           ['CONFIG_CHARM_ATSYNC']
define_memory =           ['CONFIG_USE_MEMORY']
define_projections =      ['CONFIG_USE_PROJECTIONS']
define_papi  =            ['CONFIG_USE_PAPI','PAPI3']

# Debugging defines

define_trace =            ['CELLO_TRACE']
define_debug =            ['CELLO_DEBUG']
define_debug_verbose =    ['CELLO_DEBUG_VERBOSE']

# Library defines

define_hdf5  =            ['H5_USE_16_API']
define_png   =            ['NO_FREETYPE']

# Temporary defines
define_new_initial =           ['NEW_INITIAL']

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
if (new_initial != 0):   defines = defines + define_new_initial
if (debug != 0):         defines = defines + define_debug
if (debug_verbose != 0): defines = defines + define_debug_verbose
if (memory != 0):        defines = defines + define_memory
if (atsync != 0):        defines = defines + define_atsync

#======================================================================
# ARCHITECTURE SETTINGS
#======================================================================

is_arch_valid = 0
sys.path.append("./config");
if   (arch == "linux-gnu"):    from linux_gnu    import *
elif (arch == "linux-jump"):   from linux_jump   import *
elif (arch == "ncsa-bw"):      from ncsa_bw      import *
elif (arch == "triton-pgi"):   from triton_pgi   import *
elif (arch == "triton-intel"): from triton_intel import *
elif (arch == "triton-gnu"):   from triton_gnu   import *
elif (arch == "triton-tau"):   from triton_tau   import *
elif (arch == "triton-jump"):  from triton_jump  import *

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
cc ['charm']  = charmc + charm_perf + ' '

#======================================================================
# UNIT TEST SETTINGS
#======================================================================

if (type == "serial"):
     serial_run   = ""
     parallel_run = ""
elif (type == "mpi"):
     serial_run   = ""
     parallel_run = "mpirun -np 8"
elif (type == "charm"):
     serial_run   = ""
     parallel_run = charm_path + "/bin/charmrun +p4 "

if (use_valgrind):
     valgrind = "valgrind --leak-check=full"
     parallel_run = parallel_run + " " + valgrind
     serial_run   = valgrind + " " + serial_run

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

cxxflags = flags_arch
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

linkflags    = flags_arch
linkflags    = linkflags + ' ' + flags_link
linkflags    = linkflags + ' ' + flags_config
if (type=="charm"):linkflags    = linkflags + ' ' + flags_link_charm

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
