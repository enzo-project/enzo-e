import os
import sys
import subprocess
import time
import socket

#======================================================================
# USER CONFIGURATION
#======================================================================


#----------------------------------------------------------------------
# Temporary setting for using new Output implementation
#----------------------------------------------------------------------

new_output = 0

#----------------------------------------------------------------------
# Temporary setting for using new PPM routines from enzo-dev
#----------------------------------------------------------------------

new_ppm = 1

#----------------------------------------------------------------------
# Maximum number of procesess per shared-memory node (can be larger than needed)
#----------------------------------------------------------------------

node_size = 64

#----------------------------------------------------------------------
# Whether to print out detailed messages with the TRACE() series of statements
#----------------------------------------------------------------------

trace = 0

#----------------------------------------------------------------------
# Whether to trace main phases
#----------------------------------------------------------------------

verbose = 0

#----------------------------------------------------------------------
# Whether to print out messages with the TRACE_CHARM() and TRACEPUP()
#  series of statements
#----------------------------------------------------------------------

trace_charm = 0

#----------------------------------------------------------------------
# Whether to enable displaying messages with the DEBUG() series of
# statements. Also writes messages to out.debug.<P> where P is the
# (physical) process rank. Still requires the "DEBUG" group to be
# enabled in Monitor (that is Monitor::is_active("DEBUG") must be true
# for any output)
#----------------------------------------------------------------------

debug = 0
debug_field = 0
debug_field_face = 0

#----------------------------------------------------------------------
# Do extra run-time checking.  Useful for debugging, but can potentially
# slow calculations down
#----------------------------------------------------------------------

check = 0

#----------------------------------------------------------------------
# Whether to periodically print all field values.  See
# src/Field/field_FieldBlock.cpp
#----------------------------------------------------------------------

debug_verbose = 0

#----------------------------------------------------------------------
# Whether to track dynamic memory statistics.  Can be useful, but can
# cause problems on some systems that also override new [] () / delete
# [] ()
#----------------------------------------------------------------------

memory = 1

#----------------------------------------------------------------------
# Set to 1 if Charm++ version is >= 6.7.0
#----------------------------------------------------------------------

new_charm = 1

#----------------------------------------------------------------------
# Enable charm++ dynamic load balancing
#----------------------------------------------------------------------

balance = 1

balancer = ['CommonLBs']

#----------------------------------------------------------------------
# Whether to compile with -pg to use gprof for performance profiling
#----------------------------------------------------------------------

use_gprof = 0

#----------------------------------------------------------------------
# Whether to compile with the Grackle chemistry and cooling library
#
# WARNING: must update grackle-related lines in src/Enzo/enzo.ci
#----------------------------------------------------------------------

use_grackle = 0

#----------------------------------------------------------------------
# Whether to run the test programs using valgrind to check for memory leaks
#----------------------------------------------------------------------

use_valgrind = 0

#----------------------------------------------------------------------
# Whether to use Cello Performance class for collecting performance
# data (currently requires global reductions, and may not be fully
# functional) (basic time data on root processor is still output)
#----------------------------------------------------------------------

use_performance = 1

#----------------------------------------------------------------------
# Whether to compile the CHARM++ version for use with the Projections
# performance tool.
#----------------------------------------------------------------------

use_projections = 0

#----------------------------------------------------------------------
# How many processors to run parallel unit tests
#----------------------------------------------------------------------

ip_charm = '4'

#----------------------------------------------------------------------
# Whether this is a Git repository
#----------------------------------------------------------------------

have_git = 1

#----------------------------------------------------------------------
# Whether this is a Mercurial repository
#----------------------------------------------------------------------

have_mercurial = 0

#----------------------------------------------------------------------
# Whether to use the jemalloc library for memory allocation
#----------------------------------------------------------------------

use_jemalloc = 0

#----------------------------------------------------------------------
# AUTO CONFIGURATION
#----------------------------------------------------------------------

# Whether the system has the PAPI performance API installed
# (config include may override use_papi)

use_papi = 0

#-----------------------------------------------------------------------
# COMMAND-LINE ARGUMENTS
#-----------------------------------------------------------------------

# scons command line (overrides CELLO_* environment variables)

arch = ARGUMENTS.get('arch','unknown')
prec = ARGUMENTS.get('prec','unknown')

# use environment variable if scons command line not provided

if (arch == 'unknown' and "CELLO_ARCH" in os.environ):
     arch = os.environ["CELLO_ARCH"]
if (prec == 'unknown' and "CELLO_PREC" in os.environ):
     prec = os.environ["CELLO_PREC"]

print 
print "    CELLO_ARCH scons arch=",arch
print "    CELLO_PREC scons prec=",prec
print 

#----------------------------------------------------------------------
# CONFIGURATION DEFINES
#----------------------------------------------------------------------

define = {}

# Temporary defines

# Global defines

define_cello =        ['CONFIG_USE_CELLO']

# Precision defines

define["single"] =    ['CONFIG_PRECISION_SINGLE']
define["double"] =    ['CONFIG_PRECISION_DOUBLE']
define_int_size  =    ['SMALL_INTS']

# Grackle defines

define_grackle   = ['CONFIG_USE_GRACKLE']
grackle_path     = 'grackle_path_not_set'

# Jemalloc defines
define_jemalloc  = ['CONFIG_USE_JEMALLOC']

# Performance defines

define_memory =       ['CONFIG_USE_MEMORY']
define_new_charm =    ['CONFIG_NEW_CHARM']
define_projections =  ['CONFIG_USE_PROJECTIONS']
define_performance =  ['CONFIG_USE_PERFORMANCE']
define_papi  =        ['CONFIG_USE_PAPI','PAPI3']

# Experimental code defines

define_new_output      = ['NEW_OUTPUT']
define_new_ppm         = ['NEW_PPM']

# Debugging defines

define_trace =        ['CELLO_TRACE']
define_verbose =      ['CELLO_VERBOSE']
define_trace_charm =  ['CELLO_TRACE_CHARM']
define_debug =        ['CELLO_DEBUG']
define_debug_field =  ['DEBUG_FIELD']
define_debug_field_face =  ['DEBUG_FIELD_FACE']
define_check =        ['CELLO_CHECK']

define_debug_verbose = ['CELLO_DEBUG_VERBOSE']

# Library defines

define_hdf5  =        []

define_png   =        ['NO_FREETYPE']

# Charm defines

define_charm =        ['CONFIG_USE_CHARM']  # used for Grackle 

# Python version defines

define_python_lt_27 = ['CONFIG_PYTHON_LT_27']

# Version control defines (Git or Mercurial)

define_have_version_control = ['CONFIG_HAVE_VERSION_CONTROL']


#======================================================================
# ARCHITECTURE SETTINGS
#======================================================================

is_arch_valid = 0

# Assume Python is new, but may be overridden in machine configuration
# files if needed.  For example, gordon and comet have Python 2.6
# installed, but subprocess.check_output() used below requires 2.7, so
# we set python_lt_27 = 1 in those configuration files to avoid
# calling check_output()

python_lt_27 = 0


sys.path.append("./config");

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
boost_inc = ''
boost_lib = ''
serial_run   = ""
serial_arg   = ""
parallel_run = ""
parallel_arg = ""
smp = 0

if   (arch == "gordon_gnu"):   from gordon_gnu   import *
elif (arch == "gordon_intel"): from gordon_intel import *
elif (arch == "gordon_pgi"):   from gordon_pgi   import *
elif (arch == "comet_gnu"):    from comet_gnu    import *
elif (arch == "linux_gnu"):    from linux_gnu    import *
elif (arch == "linux_intel"):  from linux_intel  import *
elif (arch == "linux_yt"):     from linux_yt     import *
elif (arch == "linux_gprof"):  from linux_gprof  import *
elif (arch == "linux_mpe"):    from linux_mpe    import *
elif (arch == "linux_tau"):    from linux_tau    import *
elif (arch == "ncsa_bw_net"):   from ncsa_bw_net import *
elif (arch == "ncsa_bw_smp"):   from ncsa_bw_smp import *
elif (arch == "faraday_gnu"):  from faraday_gnu  import *
elif (arch == "faraday_gnu_debug"):  from faraday_gnu_debug  import *
elif (arch == "mf_gnu"):       from mf_gnu       import *
elif (arch == "mf_gnu_debug"): from mf_gnu_debug import *
elif (arch == "stampede_gnu"): from stampede_gnu import *
elif (arch == "stampede_intel"): from stampede_intel import *
elif (arch == "davros_gnu"):   from davros_gnu   import *
elif (arch == "davros_gnu_debug"):  from davros_gnu_debug  import *
elif (arch == "darwin_gnu"):   from darwin_gnu   import *
elif (arch == "darwin_homebrew"):   from darwin_homebrew   import *

#======================================================================
# END ARCHITECTURE SETTINGS
#======================================================================

if (not is_arch_valid):
   print "Unrecognized architecture ",arch
   sys.exit(1)

#----------------------------------------------------------------------
# ASSEMBLE DEFINES
#----------------------------------------------------------------------

defines     = []

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

defines = defines + define_int_size

defines = defines + [{'CONFIG_NODE_SIZE' : node_size }]
defines = defines + [{'CONFIG_NODE_SIZE_3' : node_size*3 }]

defines = defines + define_hdf5
defines = defines + define_png

charm_perf = ''

if (use_projections == 1):
     defines = defines + define_projections
     charm_perf = '-tracemode projections'

if (use_performance == 1):
     defines = defines + define_performance

if (use_gprof == 1):
     flags_config = flags_config + ' -pg'

if (use_jemalloc == 1):
   defines = defines + define_jemalloc

if (use_papi != 0):      defines = defines + define_papi
if (use_grackle != 0):   defines = defines + define_grackle

if (new_output != 0):    defines = defines + define_new_output
if (new_ppm != 0):       defines = defines + define_new_ppm

if (trace != 0):         defines = defines + define_trace
if (verbose != 0):       defines = defines + define_verbose
if (trace_charm != 0):   defines = defines + define_trace_charm
if (debug != 0):         defines = defines + define_debug
if (debug_field != 0):   defines = defines + define_debug_field
if (debug_field_face != 0): defines = defines + define_debug_field_face
if (check != 0):         defines = defines + define_check
if (debug_verbose != 0): defines = defines + define_debug_verbose
if (memory != 0):        defines = defines + define_memory
if (new_charm != 0):     defines = defines + define_new_charm
if (python_lt_27 != 0):  defines = defines + define_python_lt_27
if (have_git != 0 or have_mercurial != 0 ):defines = defines + define_have_version_control

defines = defines + define_charm
defines = defines + define_cello

#======================================================================
# FINAL CHARM SETUP
#======================================================================

charmc = charm_path + '/bin/charmc -language charm++ '

cxx = charmc + charm_perf + ' '

if (balance == 1):
     flags_cxx_charm = flags_cxx_charm + " -balancer " + " -balancer ".join(balancer)
     flags_link_charm = flags_link_charm + " -module " + " -module ".join(balancer)

#======================================================================
# UNIT TEST SETTINGS
#======================================================================

if (parallel_run == ''):
   if (smp == 1):
      parallel_run = charm_path + "/bin/charmrun ++ppn " + ip_charm + " +p" + ip_charm
   else:
      parallel_run = charm_path + "/bin/charmrun +p" + ip_charm

if (use_valgrind):
     valgrind = "valgrind --leak-check=full"
     serial_run   = valgrind + " " + serial_run
     parallel_run = parallel_run + " " + valgrind

#======================================================================
# CELLO PATHS
#======================================================================

bin_path = '#/bin'
lib_path = '#/lib'
inc_path = '#/include'
test_path= 'test'

Export('bin_path')
Export('grackle_path')
Export('use_grackle')
Export('use_jemalloc')
Export('lib_path')
Export('inc_path')
Export('test_path')
Export('ip_charm')
Export('smp')


cpppath     = [inc_path]
fortranpath = [inc_path]
libpath     = [lib_path]

#----------------------------------------------------------------------
# PAPI PATHS
#----------------------------------------------------------------------

if (use_papi):
     cpppath = cpppath + [papi_inc]
     libpath = libpath + [papi_lib]

#----------------------------------------------------------------------
# HDF5 PATHS
#----------------------------------------------------------------------

cpppath = cpppath + [ hdf5_inc ]
libpath = libpath + [ hdf5_lib ]

#----------------------------------------------------------------------
# BOOST PATHS
#----------------------------------------------------------------------

cpppath = cpppath + [ boost_inc ]
libpath = libpath + [ boost_lib ]

#----------------------------------------------------------------------
# GRACKLE PATH
#----------------------------------------------------------------------

if (use_grackle != 0):
      cpppath.append(grackle_path)
      libpath.append(grackle_path)

#----------------------------------------------------------------------
# LIBPNG PATHS
#----------------------------------------------------------------------

cpppath = cpppath + [ png_path + '/include' ]
libpath = libpath + [ png_path + '/lib']

#----------------------------------------------------------------------
# FORTRAN LINK PATH
#----------------------------------------------------------------------

libpath = libpath + [libpath_fortran]


#======================================================================
# ENVIRONMENT
#======================================================================

environ  = os.environ

cxxflags = flags_arch + ' ' + flags_arch_cpp
cxxflags = cxxflags + ' ' + flags_cxx
cxxflags = cxxflags + ' ' + flags_config
cxxflags = cxxflags + ' ' + flags_cxx_charm
Export('cxxflags')

cflags   = flags_arch
cflags   = cflags + ' ' + flags_cc
cflags   = cflags + ' ' + flags_config
cflags   = cflags + ' ' + flags_cc_charm

fortranflags = flags_arch
fortranflags = fortranflags + ' ' + flags_fc
fortranflags = fortranflags + ' ' + flags_config
fortranflags = fortranflags + ' ' + flags_fc_charm

linkflags    = flags_arch + ' ' + flags_arch_cpp
linkflags    = linkflags + ' ' + flags_link
linkflags    = linkflags + ' ' + flags_config
linkflags    = linkflags + ' ' + flags_link_charm

if (prec == 'double'):
    fortranflags = fortranflags + ' ' + flags_prec_double
if (prec == 'single'):
    fortranflags = fortranflags + ' ' + flags_prec_single

if not os.path.exists("include"):
     os.makedirs("include")
cello_def = open ("include/auto_config.def", "w")


env = Environment (
     CC           = cc,
     CFLAGS       = cflags,
     CPPDEFINES   = defines,
     CPPPATH      = cpppath,
     CXX          = cxx,
     CXXFLAGS     = cxxflags,
     ENV          = environ,
     FORTRANFLAGS = fortranflags,
     FORTRAN      = f90,
     FORTRANLIBS  = libs_fortran,
     FORTRANPATH  = fortranpath,
     LIBPATH      = libpath,
     LINKFLAGS    = linkflags )

cello_def.write ("#define CELLO_ARCH "
		"\""+arch+"\"\n")
cello_def.write ("#define CELLO_PREC "
		"\""+prec+"\"\n")
cello_def.write ("#define CELLO_CC "
		"\""+cc+"\"\n")	
cello_def.write ("#define CELLO_CFLAGS "
		"\""+cflags+"\"\n")
cello_def.write ("#define CELLO_CPPDEFINES "
		"\""+" ".join(map(str,defines))+"\"\n")
cello_def.write ("#define CELLO_CPPPATH "
		"\""+" ".join(map(str,cpppath))+"\"\n")
cello_def.write ("#define CELLO_CXX "
		"\""+cxx+"\"\n")	
cello_def.write ("#define CELLO_CXXFLAGS "
		"\""+cxxflags+"\"\n")
cello_def.write ("#define CELLO_FORTRANFLAGS "
		"\""+fortranflags+"\"\n")
cello_def.write ("#define CELLO_FORTRAN "
		"\""+f90+"\"\n")
cello_def.write ("#define CELLO_FORTRANLIBS "
		"\""+" ".join(map(str,libs_fortran))+"\"\n")
cello_def.write ("#define CELLO_FORTRANPATH "
		"\""+" ".join(map(str,fortranpath))+"\"\n")
cello_def.write ("#define CELLO_LIBPATH "
		"\""+" ".join(map(str,libpath))+"\"\n")
cello_def.write ("#define CELLO_LINKFLAGS "
		"\""+linkflags+"\"\n" )
cello_def.write ("#define CELLO_HOST "
		"\""+socket.gethostname()+"\"\n" )
cello_def.write ("#define CELLO_DIR "
		"\""+Dir('.').abspath+"\"\n" )
cello_def.write ("#define CELLO_DATE "
		"\""+time.strftime("%Y-%m-%d",time.gmtime())+"\"\n" )
cello_def.write ("#define CELLO_TIME "
		"\""+time.strftime("%H:%M:%S",time.gmtime())+"\"\n" )

#----------
# Python version >= 2.7 is required for subprocess.check_output()

if (python_lt_27 == 0):
     charm_version =  subprocess.check_output (["cat", charm_path + "/VERSION"]).rstrip();
     cello_def.write ("#define CELLO_CHARM_VERSION "+charm_version+"\n" )
     
     fp_charm_version = open ("test/CHARM_VERSION", "w")
     fp_charm_version.write(charm_version + "\n");
     fp_charm_version.close()
     		      
else:
     cello_def.write ("#define CELLO_CHARM_VERSION 0\n")	
     fp_charm_version = open ("test/CHARM_VERSION", "w")
     fp_charm_version.write("unknown\n");
     fp_charm_version.close()

Clean('.','test/CHARM_VERSION')

cello_def.write ("#define CELLO_CHARM_PATH \"" + charm_path + "\"\n" )

#----------
# Both Python version 2.7 is required, and git must be installed

if (python_lt_27 == 0 and have_git):

     cello_def.write ("#define CELLO_CHANGESET "
		      "\""+subprocess.check_output
		      (["git", "rev-parse", "HEAD"]).rstrip()+"\"\n" )

elif (python_lt_27 == 0 and have_mercurial):

     cello_def.write ("#define CELLO_CHANGESET "
                      "\""+subprocess.check_output
                      (["hg","id","-n"]).rstrip()+"\"\n" )

else:
     cello_def.write ("#define CELLO_CHANGESET \"unknown\"\n" )

#----------

# Find how Charm++ was compiled using the bin symbolic link real path

t = os.path.realpath(charm_path + '/bin')
i0=t.rfind('/')
i0=t.rfind('/',0,i0)
i1=t.rfind('/bin')
charm_build = t[i0+1:i1]

cello_def.write ("#define CHARM_BUILD \"" + charm_build + "\"\n")
fp_charm_build = open ("test/CHARM_BUILD", "w")
fp_charm_build.write(charm_build + "\n");
fp_charm_build.close()
cello_def.close()
Clean('.','test/CHARM_BUILD')

#======================================================================
# BUILDERS
#======================================================================

charm_builder = Builder (action="${CXX} $SOURCE; mv ${ARG}.*.h `dirname $SOURCE`")
cpp_builder = Builder (action="/usr/bin/cpp -E $_CPPDEFFLAGS $SOURCE > $TARGET")
env.Append(BUILDERS = { 'CharmBuilder' : charm_builder })
env.Append(BUILDERS = { 'CppBuilder'   : cpp_builder })

Export('env')
Export('parallel_run')
Export('parallel_arg')
Export('serial_run')
Export('serial_arg')
Export('use_papi')

# Build in build-<branch> directory if this is a git repository

if (have_git == 1):
   branch = subprocess.check_output(['git', 'rev-parse', '--abbrev-ref', 'HEAD']).rstrip()
   build_dir = 'build-' + branch
else:     
   build_dir = 'build'
   
SConscript( 'src/SConscript',variant_dir=build_dir)
SConscript('test/SConscript')

#======================================================================
# CLEANING
#======================================================================

# non-permanent directories
Clean('.','bin')
Clean('.','lib')
Clean('.','src-html')
Clean('.','src-latex')
Clean('.','src-xml')

# files left behind by enzo-p
Clean('.','Checkpoint')
Clean('.','parameters.out')
Clean('.','parameters.libconfig')

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
