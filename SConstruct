import os
import sys

#----------------------------------------------------------------------
# USER CONFIGURATION
#----------------------------------------------------------------------

# Whether to print out messages with the TRACE() statements

trace           = 1

# Whether to periodically print basic statistics about field values.
# See src/Field/field_FieldBlock.cpp

debug           = 0

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

# Temporary code variation for testing purposes.  Changing the value
# below may result in the code hanging, crashing, or producing
# incorrect results.

original_refresh = 0

# Triton MPI type (openmpi or mpich2)

mpi_type = 'mpich2'

#----------------------------------------------------------------------
# AUTO CONFIGURATION
#----------------------------------------------------------------------

# Whether the system has the PAPI performance API installed

use_papi = 0

env = Environment()

if not env.GetOption('clean'):

     conf = Configure(env)

     if not conf.CheckCHeader('papi.h'):
          print 'PAPI not installed'
          use_papi = 0
     else:
          print 'PAPI installed'
          use_papi = 1

     env = conf.Finish()

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

# Parallel type defines
define = {}

define["serial"] =        []
define["mpi"]    =        ['CONFIG_USE_MPI']
define["charm"]  =        ['CONFIG_USE_CHARM']

# Precision defines
define["single"] =        ['CONFIG_PRECISION_SINGLE']
define["double"] =        ['CONFIG_PRECISION_DOUBLE']

# Library defines
define_hdf5  =            ['H5_USE_16_API']
define_png   =            ['NO_FREETYPE']
define_papi  =            ['CONFIG_USE_PAPI']

# Debugging defines
define_trace =            ['CELLO_TRACE']
define_debug =            ['CELLO_DEBUG']
define_debug_verbose =    ['CELLO_DEBUG_VERBOSE']

# Performance defines
define_atsync =           ['CONFIG_CHARM_ATSYNC']
define_memory =           ['CONFIG_USE_MEMORY']
define_projections =      ['CONFIG_USE_PROJECTIONS']

# Experimental code defines
define_original_refresh = ['ORIGINAL_REFRESH']

defines     = []
defines_xlc = ""
defines_xlf = ""

#----------------------------------------------------------------------
# ASSEMBLE DEFINES
#----------------------------------------------------------------------

# IBM's inconsistent -D defines for C and Fortran can break scons: DC
# and DF are workarounds

DC = ' -D'
DF = ' -WF,-D'

# Parallel type

if (type == 'serial'):
     defines = defines
elif (type == 'mpi' or type == 'charm'):
     defines     = defines          + define[type]
     defines_xlc = defines_xlc + DC + define[type][0]
     defines_xlf = defines_xlf + DF + define[type][0]
else:
     print "Unrecognized parallel type ",type
     print
     print "Valid types are 'serial', 'mpi', and 'charm'"
     print
     print "The type is set using the environment variable $CELLO_TYPE"
     print "or by using 'scons type=<type>"
     sys.exit(1)

# Precision

if (prec == 'single' or prec == 'double'):
     defines     = defines                 + define[prec]
     defines_xlc = defines_xlc + DC     + define[prec][0]
     defines_xlf = defines_xlf + DF + define[prec][0]
else:
     print "Unrecognized precision ",prec
     print
     print "Valid precisions are 'single' and 'double'"
     print
     print "The precision is set using the environment variable $CELLO_PREC"
     print "or by using 'scons prec=<precision>"
     sys.exit(1)

# CHARM++ Projections

charm_perf = ''

if (use_projections == 1):
     charm_perf = '-tracemode projections'
     defines     = defines          + define_projections
     defines_xlc = defines_xlc + DC + define_projections[0]
     defines_xlf = defines_xlf + DF + define_projections[0]

# GProf flags

flags_gprof = ''

if (use_gprof == 1):
     flags_gprof = '-pg '
     
# Experimental defines

if (original_refresh == 1):
     defines    =  defines          + define_original_refresh
     defines_xlc = defines_xlc + DC + define_original_refresh[0]
     defines_xlf = defines_xlf + DF + define_original_refresh[0]

# PAPI defines

if (use_papi != 0): 
     defines     = defines          + define_papi
     defines_xlc = defines_xlc + DC + define_papi[0]
     defines_xlf = defines_xlf + DF + define_papi[0]

# TRACE defines

if (trace != 0):
     defines     = defines          + define_trace
     defines_xlc = defines_xlc + DC + define_trace[0]
     defines_xlf = defines_xlf + DF + define_trace[0]

# Debug defines

if (debug != 0):
     defines     = defines          + define_debug
     defines_xlc = defines_xlc + DC + define_debug[0]
     defines_xlf = defines_xlf + DF + define_debug[0]

# debug_verbose defines

if (debug_verbose != 0):
     defines     = defines          + define_debug_verbose
     defines_xlc = defines_xlc + DC + define_debug_verbose[0]
     defines_xlf = defines_xlf + DF + define_debug_verbose[0]

# memory defines

if (memory != 0):
     defines     = defines          + define_memory
     defines_xlc = defines_xlc + DC + define_memory[0]
     defines_xlf = defines_xlf + DF + define_memory[0]

# atsync defines

if (atsync != 0):
     defines     = defines          + define_atsync
     defines_xlc = defines_xlc + DC + define_atsync[0]
     defines_xlf = defines_xlf + DF + define_atsync[0]

# HDF5 library defines

defines     = defines     +         define_hdf5
defines_xlc = defines_xlc + DC + define_hdf5[0] 
defines_xlf = defines_xlf + DF + define_hdf5[0] 

# PNG library defines

defines     = defines             + define_png
defines_xlc = defines_xlc + DC + define_png[0]
defines_xlf = defines_xlf + DF + define_png[0]

#======================================================================
# ARCHITECTURE SETTINGS
#======================================================================

fortran = {}
cxx     = {}
cc      = {}

cflags_arch       = ''
cxxflags_arch     = ''
fortranflags_arch = ''

path_list = []

is_arch_valid = 0

#----------------------------------------------------------------------
if (arch == "linux-gnu"):
#----------------------------------------------------------------------

     is_arch_valid = 1

     charm_path = '/home/bordner/charm/charm'  # arch

     fortran['serial'] = 'gfortran'
     fortran['mpi']    = 'gfortran'
     fortran['charm']  = 'gfortran'

     cxx['serial'] = 'g++'
     cxx['mpi']    = 'mpic++'

     cc['serial']  = 'gcc'
     cc['mpi']     = 'mpicc'

# for architecture-dependent defines

     cppdefines        = defines

# for extra Fortran libraries

     libpath_fortran = ''
     libs_fortran    = ['gfortran']
     # -rdynamic is to include symbols in backtrace
     fortranlinkflags_arch  = '-rdynamic'

# PAPI path (optional)
     papi_path = '/usr/local'
# HDF5 path
     hdf5_path = '/usr'
# Optional debugging flags
     flags_debug = '-g'
# Optional optimization flags
     flags_opt   = ''
# Optional warnings-level flags
     flags_warn  = '-Wall'

#----------------------------------------------------------------------
elif (arch == "ncsa-bd"):
#----------------------------------------------------------------------

     is_arch_valid = 1

     charm_path = '/home/bordner/charm/charm'

     fc_path  = '/opt/ibmcmp/xlf/13.1'
     cc_path  = '/opt/ibmcmp/vac/11.1'
     cxx_path = '/opt/ibmcmp/vacpp/11.1'

     fortran['serial'] = fc_path + '/bin/xlf_r'
     fortran['mpi']    = fc_path + '/bin/xlf_r'
     fortran['charm']  = fc_path + '/bin/xlf_r'

     cxx['serial'] = cxx_path + '/bin/xlC_r'
     cxx['mpi']    = 'mpCC'

     cc['serial']  = cc_path + '/bin/xlc_r'
     cc['mpi']     = 'mpcc'

# Architecture-dependent flags

     cflags_arch       = defines_xlc
     cxxflags_arch     = defines_xlc
     fortranflags_arch = defines_xlf + '-qextname'

# Extra fortran libraries

     libpath_fortran = fc_path + '/lib64'
     libs_fortran    = ['xlf90','xlfmath','xl']
     fortranlinkflags_arch  = ''

# PAPI path (optional)
     papi_path = '/opt/usersoft/papi/4.1.0'
# HDF5 path
     hdf5_path = '/home/bordner'
# Optional debugging flags
     flags_debug = ''
# Optional optimization flags
     flags_opt   = '-O3 -qhot -q64'
# Optional warnings-level flags
     flags_warn  = ''

#----------------------------------------------------------------------
elif (arch == "triton-pgi"):
#----------------------------------------------------------------------

     is_arch_valid = 1

     # Requires modules pgi, mpich_mx

     charm_path = '/home/jobordner/public/charm/charm-' + mpi_type + '-pgi'

     fortran['serial'] = 'pgf90'
     fortran['mpi']    = 'pgf90'
     fortran['charm']  = 'pgf90'

     cxx['serial'] = 'pgCC'
     cxx['mpi']    = 'mpicxx'

     cc['serial']  = 'pgcc'
     cc['mpi']     = 'mpicc'

# Architecture-dependent defines

     cppdefines        = defines

# For extra fortran libraries

     libpath_fortran = ''
     libs_fortran    = []
     # -rdynamic is to include symbols in backtrace
     fortranlinkflags_arch  = '-pgf90libs'

# PAPI path (optional)
     papi_path = ''
# HDF5 path
     hdf5_path = '/opt/hdf5/pgi'
# Optional debugging flags
     flags_debug = '-Ktrap=fp -g'
# Optional optimization flags
     flags_opt   = ''
# Optional warnings-level flags
     flags_warn  = ''

#----------------------------------------------------------------------
elif (arch == "triton-intel"):
#----------------------------------------------------------------------

     is_arch_valid = 1

     # Requires modules intel mpich_mx

     charm_path = '/home/jobordner/public/charm/charm-' + mpi_type + '-intel'

     fortran['serial'] = 'ifort'
     fortran['mpi']    = 'ifort'
     fortran['charm']  = 'ifort'

     cxx['serial'] = 'icpc'
     cxx['mpi']    = 'mpicxx'

     cc['serial']  = 'icc'
     cc['mpi']     = 'mpicc'

# Architecture-dependent defines

     cppdefines        = defines

# For extra fortran libraries

     libpath_fortran = ''
     libs_fortran    = ['imf','ifcore','ifport','stdc++']
     fortranlinkflags_arch  = ''

# PAPI path (optional)
     papi_path = ''
# HDF5 path
     hdf5_path = '/opt/hdf5/intel'
# Optional debugging flags
     flags_debug = '-g'
# Optional optimization flags
     flags_opt   = ''
# Optional warnings-level flags
     flags_warn  = ''

#----------------------------------------------------------------------
elif (arch == "triton-gnu"):
#----------------------------------------------------------------------

     is_arch_valid = 1

     # Requires modules gnu, mpich_mx

     charm_path = '/home/jobordner/public/charm/charm-' + mpi_type + '-gnu'

     fortran['serial'] = 'gfortran'
     fortran['mpi']    = 'gfortran'
     fortran['charm']  = 'gfortran'

     cxx['serial'] = 'g++'
     cxx['mpi']    = 'mpic++'

     cc['serial']  = 'gcc'
     cc['mpi']     = 'mpicc'

# Architecture-dependent defines

     cppdefines        = defines

# For extra fortran libraries

     libpath_fortran = ''
     libs_fortran    = ['gfortran']
     fortranlinkflags_arch  = '-rdynamic'

# PAPI path (optional)
     papi_path = ''
# HDF5 path
     hdf5_path = '/opt/hdf5/gnu'
# Optional debugging flags
     flags_debug = '-g'
# Optional optimization flags
#     flags_opt   = '-O3'
     flags_opt   = ''
# Optional warnings-level flags
     flags_warn  = '-Wall'

#======================================================================
# END ARCHITECTURE SETTINGS
#======================================================================

if (not is_arch_valid):
   print "Unrecognized architecture type ",arch
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
     parallel_run = charm_path + "/bin/charmrun +p8 "

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

cxxflags = \
    cxxflags_arch + ' ' + \
    flags_debug + ' ' + \
    flags_opt   + ' ' + \
    flags_warn  + ' ' + \
    flags_gprof

cflags = \
    cflags_arch + ' ' + \
    flags_debug + ' ' + \
    flags_opt   + ' ' + \
    flags_gprof  + \
    flags_warn

fortranflags = \
    fortranflags_arch + ' ' + \
    flags_debug + ' ' + \
    flags_opt   + ' ' + \
    flags_warn  + ' ' + \
    flags_gprof

linkflags = \
    fortranlinkflags_arch  + ' ' + \
    flags_debug + ' ' + \
    flags_opt   + ' ' + \
    flags_gprof     + \
    flags_warn

env = Environment (
     CC           = cc[type],	
     CFLAGS       = cflags,
     CPPDEFINES   = cppdefines,
     CPPPATH      = cpppath,
     CXX          = cxx[type],	
     CXXFLAGS     = cxxflags,
     ENV          = environ,
     FORTRANFLAGS = fortranflags,
     FORTRAN      = fortran[type],
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
