import os
import sys

#----------------------------------------------------------------------
# CONFIGURATION
#----------------------------------------------------------------------

trace           = 0
debug           = 0
debug_verbose   = 0

atsync          = 0 # CHARM++: only load balance when AtSync() called
trace           = 0

use_gprof       = 0
use_papi        = 0
use_valgrind    = 0
use_projections = 1

# temporary code variations for testing
original_refresh = 1
skip_reduce      = 0  # constant dt, which is currently hard-coded!

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

define_hdf5  =  ['H5_USE_16_API'];
define_png   =  ['NO_FREETYPE'];
define_papi  =  ['CONFIG_USE_PAPI'];
define_trace =  ['CELLO_TRACE'];
define_atsync =  ['CONFIG_CHARM_ATSYNC'];
define_debug =  ['CELLO_DEBUG'];
define_debug_verbose =  ['CELLO_DEBUG_VERBOSE'];
define_projections =  ['CONFIG_USE_PROJECTIONS']
define_original_refresh = ['ORIGINAL_REFRESH']

define_skip_reduce = ['TEMP_SKIP_REDUCE']

defines     = []
defines_xlc = ""
defines_xlf = ""

#--------------------------------------------------
# PARALLEL DEFINES
#--------------------------------------------------

if (type == 'serial'):
     defines = defines
elif (type == 'mpi'):
     defines     = defines              + define_mpi
     defines_xlc = defines_xlc + ' -D'  + define_mpi[0]
     defines_xlf = defines_xlf + ' -WF,-D'  + define_mpi[0]
elif (type == 'charm'):
     defines     = defines              + define_charm
     defines_xlc = defines_xlc + ' -D'  + define_charm[0]
     defines_xlf = defines_xlf + ' -WF,-D'  + define_charm[0]
else:
     print "Unrecognized parallel type ",type
     print
     print "Valid types are 'serial', 'mpi', and 'charm'"
     print
     print "The type is set using the environment variable $CELLO_TYPE"
     print "or by using 'scons type=<type>"
     sys.exit(1)

#======================================================================
# CHARM++ PROJECTIONS
#======================================================================


charm_perf        = ''
cxxflags_perf     = ''
cflags_perf       = ''
fortranflags_perf = ''
linkflags_perf    = ''

if (use_projections == 1):
     charm_perf = '-tracemode projections'
     defines     = defines              + define_projections
     defines_xlc = defines_xlc + ' -D'  + define_projections[0]
     defines_xlf = defines_xlf + ' -WF,-D'  + define_projections[0]

flags_gprof = ''

if (use_gprof == 1):
     flags_gprof = '-pg '
     
#--------------------------------------------------

if (prec == 'single'):
     defines = defines + define_single
     defines_xlc = defines_xlc + ' -D' + define_single[0]
     defines_xlf = defines_xlf + ' -WF,-D' + define_single[0]
elif (prec == 'double'):
     defines = defines + define_double
     defines_xlc = defines_xlc + ' -D' + define_double[0]
     defines_xlf = defines_xlf + ' -WF,-D' + define_double[0]
else:
     print "Unrecognized precision ",prec
     print
     print "Valid precisions are 'single' and 'double'"
     print
     print "The precision is set using the environment variable $CELLO_PREC"
     print "or by using 'scons prec=<precision>"
     sys.exit(1)

print defines
print defines_xlc
print defines_xlf

#-----------------------------------------------------------------------
if (original_refresh == 1):
     defines = defines + define_original_refresh
     defines_xlc = defines_xlc + ' -D' + define_original_refresh[0]
     defines_xlf = defines_xlf + ' -WF,-D' + define_original_refresh[0]

if (skip_reduce == 1):
     defines = defines + define_skip_reduce
     defines_xlc = defines_xlc + ' -D' + define_skip_reduce[0]
     defines_xlf = defines_xlf + ' -WF,-D' + define_skip_reduce[0]

#-----------------------------------------------------------------------
# Display configuration settings
#-----------------------------------------------------------------------

print "CONFIGURATION"
print 
print "    CELLO_ARCH: scons arch=",arch
print "    CELLO_TYPE: scons type=",type
print "    CELLO_PREC: scons prec=",prec
print 

#==================================================
# Initialize environment according to configuration
#==================================================
#--------------------------------------------------

if (use_papi != 0): 
     defines     = defines             + define_papi
     defines_xlc = defines_xlc + ' -D' + define_papi[0]
     defines_xlf = defines_xlf + ' -WF,-D' + define_papi[0]

if (trace != 0):
     defines = defines + define_trace
     defines_xlc = defines_xlc + ' -D' + define_trace[0]
     defines_xlf = defines_xlf + ' -WF,-D' + define_trace[0]

if (debug != 0):
     defines = defines + define_debug
     defines_xlc = defines_xlc + ' -D' + define_debug[0]
     defines_xlf = defines_xlf + ' -WF,-D' + define_debug[0]

if (debug_verbose != 0):
     defines = defines + define_debug_verbose
     defines_xlc = defines_xlc + ' -D' + define_debug_verbose[0]
     defines_xlf = defines_xlf + ' -WF,-D' + define_debug_verbose[0]

if (atsync != 0):
     defines = defines + define_atsync
     defines_xlc = defines_xlc + ' -D' + define_atsync[0]
     defines_xlf = defines_xlf + ' -WF,-D' + define_atsync[0]

#--------------------------------------------------

defines     = defines     +         define_hdf5
defines_xlc = defines_xlc + ' -D' + define_hdf5[0] 
defines_xlf = defines_xlf + ' -WF,-D' + define_hdf5[0] 

#--------------------------------------------------

defines     = defines             + define_png;
defines_xlc = defines_xlc + ' -D' + define_png[0]
defines_xlf = defines_xlf + ' -WF,-D' + define_png[0]

#======================================================================
# ARCHITECTURE SETTINGS
#======================================================================

#----------------------------------------------------------------------
if (arch == "linux"):
#----------------------------------------------------------------------

     charm_path = '/home/bordner/charm/charm'  # arch

     fortran_serial = 'gfortran'
     fortran_mpi    = 'gfortran'
     fortran_charm  = 'gfortran'

     cxx_serial = 'g++'
     cxx_mpi    = 'mpic++'
     cxx_charm  = charm_path + '/bin/charmc -language charm++ ' + charm_perf + ' '

     cc_serial  = 'gcc'
     cc_mpi     = 'mpicc'
     cc_charm   = charm_path + '/bin/charmc -language charm++ ' + charm_perf + ' '

     cppdefines = defines
     cxxflags_define     = ''
     fortranflags_define = ''
     fortranpath_lib = ''
     fortranlibs = ['gfortran']

     papi_path = '/usr/local'
     papi_inc = (papi_path + '/include')
     papi_lib = (papi_path + '/lib')


     hdf5_path = '/usr'
     hdf5_inc = (hdf5_path + '/include')
     hdf5_lib = (hdf5_path + '/lib')

     flags_debug = '-g'
#     flags_opt   = '-O3'
     flags_opt   = ''
     flags_prec  = '-m128bit-long-double'
     flags_warn  = '-Wall'

     cxxflags_debug = flags_debug
     cxxflags_opt   = flags_opt
     cxxflags_prec  = flags_prec
     cxxflags_warn  = flags_warn

     cflags_debug = flags_debug
     cflags_opt   = flags_opt
     cflags_prec  = flags_prec
     cflags_warn  = flags_warn

     fortranflags_debug = flags_debug
     fortranflags_opt   = flags_opt
     fortranflags_prec  = flags_prec
     fortranflags_warn  = flags_warn

     linkflags_arch = ''
     linkflags_debug = flags_debug
     linkflags_opt   = flags_opt
     linkflags_prec  = flags_prec
     linkflags_warn  = flags_warn

#----------------------------------------------------------------------
elif (arch == "ncsa-bd"):
#----------------------------------------------------------------------

     charm_path = '/home/bordner/charm/charm'

     fc_path  = '/opt/ibmcmp/xlf/13.1'
     cc_path  = '/opt/ibmcmp/vac/11.1'
     cxx_path = '/opt/ibmcmp/vacpp/11.1'

     fortran_serial = fc_path + '/bin/xlf_r'
     fortran_mpi    = fc_path + '/bin/xlf_r'
     fortran_charm  = fc_path + '/bin/xlf_r'

     cxx_serial = cxx_path + '/bin/xlC_r'
     cxx_mpi    = 'mpCC'
     cxx_charm  = charm_path + '/bin/charmc -language charm++ ' + charm_perf + ' '

     cc_serial  = cc_path + '/bin/xlc_r'
     cc_mpi     = 'mpcc'
     cc_charm   = charm_path + '/bin/charmc -language charm++ ' + charm_perf + ' '

# defines moved to flags since xlf_r expects -WF,-Dblah but xlC expects -Dblah

     cppdefines = ''
     cxxflags_define     = defines_xlc
     fortranflags_define = defines_xlf
     fortranpath_lib = fc_path + '/lib64'
     fortranlibs = ['xlf90','xlfmath','xl']

     papi_path = '/opt/usersoft/papi/4.1.0'
     papi_inc = (papi_path + '/include')
     papi_lib = (papi_path + '/lib64')

     hdf5_path = '/opt/hdf5-1.8.4-patch1-64bit'
     hdf5_inc = (hdf5_path + '/include')
     hdf5_lib = (hdf5_path + '/lib')

     flags_debug = ''
     flags_opt   = '-O3 -qhot -q64'
     flags_prec  = ''
     flags_warn  = ''

     cxxflags_debug = flags_debug
     cxxflags_opt   = flags_opt
     cxxflags_prec  = flags_prec
     cxxflags_warn  = flags_warn

     cflags_debug = flags_debug
     cflags_opt   = flags_opt
     cflags_prec  = flags_prec
     cflags_warn  = flags_warn

     fortranflags_debug = flags_debug + ' -qextname'
     fortranflags_opt   = flags_opt
     fortranflags_prec  = flags_prec
     fortranflags_warn  = flags_warn

     linkflags_arch = ''
     linkflags_debug = flags_debug
     linkflags_opt   = flags_opt
     linkflags_prec  = flags_prec
     linkflags_warn  = flags_warn

#----------------------------------------------------------------------
elif (arch == "sdsc-triton"):
#----------------------------------------------------------------------

     charm_path = '/home/jobordner/public/charm/charm'

     fc_path  = '/opt/openmpi/pgi/mx'
     cc_path  = '/opt/openmpi/pgi/mx'
     cxx_path = '/opt/openmpi/pgi/mx'

     fortran_serial = 'pgf90'
     fortran_mpi    = 'pgf90'
     fortran_charm  = 'pgf90'

     cxx_serial = 'pgCC'
     cxx_mpi    = 'mpicxx'
     cxx_charm  = charm_path + '/bin/charmc -language charm++ ' + charm_perf + ' '

     cc_serial  = 'pgcc'
     cc_mpi     = 'mpicc'
     cc_charm  = charm_path + '/bin/charmc -language charm++ ' + charm_perf + ' '

     cppdefines = defines
     cxxflags_define     = ''
     fortranflags_define = ''
     fortranpath_lib = ''
     fortranlibs = []

     papi_path = ''
     papi_inc = (papi_path + '/include')
     papi_lib = (papi_path + '/lib')

     hdf5_path = '/opt/hdf5/pgi'
     hdf5_inc = (hdf5_path + '/include')
     hdf5_lib = (hdf5_path + '/lib')

     flags_debug = ''
     flags_opt   = '-fast'
     flags_prec  = ''
     flags_warn  = ''

     cxxflags_debug = flags_debug
     cxxflags_opt   = flags_opt
     cxxflags_prec  = flags_prec
     cxxflags_warn  = flags_warn

     cflags_debug = flags_debug
     cflags_opt   = flags_opt
     cflags_prec  = flags_prec
     cflags_warn  = flags_warn

     fortranflags_debug = flags_debug
     fortranflags_opt   = flags_opt
     fortranflags_prec  = flags_prec
     fortranflags_warn  = flags_warn

     linkflags_arch = '-pgf90libs'
     linkflags_debug = flags_debug
     linkflags_opt   = flags_opt
     linkflags_prec  = flags_prec
     linkflags_warn  = flags_warn



#======================================================================
# PARALLELISM SETTINGS
#======================================================================

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

#======================================================================
# CELLO PATHS
#======================================================================

cpppath     = ['#/include'];
fortranpath = ['#/include'];
libpath     = ['#/lib'];

#----------------------------------------------------------------------
# PAPI PATHS
#----------------------------------------------------------------------

if (use_papi):
     cpppath = cpppath + [papi_inc]
     libpath = libpath + [papi_lib]

#----------------------------------------------------------------------
# HDF5 PATHS
#----------------------------------------------------------------------

cpppath = cpppath + [hdf5_inc]
libpath = libpath + [hdf5_lib]

#----------------------------------------------------------------------
# FORTRAN LINK PATH
#----------------------------------------------------------------------

libpath = libpath + [fortranpath_lib]

#======================================================================
# VALGRIND SETTINGS
#======================================================================

if (use_valgrind):
     valgrind = "valgrind --leakcheck=full"
     parallel_run = parallel_run + " " + valgrind
     serial_run   = vagrind + " " + serial_run

#======================================================================
# ENVIRONMENT
#======================================================================

environ  = os.environ

cxxflags = \
    cxxflags_debug + ' ' + \
    cxxflags_opt   + ' ' + \
    cxxflags_prec  + ' ' + \
    cxxflags_perf  + ' ' + \
    cxxflags_warn  + ' ' + \
    flags_gprof    + \
    cxxflags_define

cflags = \
    cflags_debug + ' ' + \
    cflags_opt   + ' ' + \
    cflags_prec  + ' ' + \
    flags_gprof  + \
    cflags_warn

fortranflags = \
    fortranflags_debug + ' ' + \
    fortranflags_opt   + ' ' + \
    fortranflags_prec  + ' ' + \
    fortranflags_warn  + ' ' + \
    flags_gprof        + \
    fortranflags_define

linkflags = \
    linkflags_arch  + ' ' + \
    linkflags_debug + ' ' + \
    linkflags_opt   + ' ' + \
    linkflags_perf  + ' ' + \
    linkflags_prec  + ' ' + \
    flags_gprof     + \
    linkflags_warn

env = Environment (
     CC           = cc,	
     CFLAGS       = cflags,
     CPPDEFINES   = cppdefines,
     CPPPATH      = cpppath,
     CXX          = cxx,	
     CXXFLAGS     = cxxflags,
     ENV          = environ,
     FORTRANFLAGS = fortranflags,
     FORTRAN      = fortran,
     FORTRANLIBS  = fortranlibs,
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

