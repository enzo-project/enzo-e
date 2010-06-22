import os
import sys

TOPDIR  = GetLaunchDir()

# Initialize platform

platform='unknown'

platform = ARGUMENTS.get('platform','unknown')

if (platform == 'unknown' and "CELLO_PLATFORM" in os.environ):
     platform = os.environ["CELLO_PLATFORM"]

#==================================================
# Check code for revision updates and local modifications
#==================================================

repository = "svn+ssh://client65-88.sdsc.edu/usr/local/svn/cello/trunk"
get_revision = "| awk '/Revision:/ {print $2}'"

# TODO: TEST FOR INTERNET CONNECTION

# revision_new     = int(os.popen("svn info " + repository + get_revision).read())
revision_current = int(os.popen("svn info " + get_revision).read())
revision_changes = int(os.popen("svn status | grep -v '?' | wc -l").read())

if (revision_changes != 0):
   print "\nWARNING: Working directory has local modifications\n"
# elif (revision_new != revision_current):
#    print "\nWARNING: Working directory is not up-to-date with the repository\n"
else:
   print "\nWorking directory synched with repository\n"

# --------------------------------------------------

#==================================================
# Initialize environment according to platform
#==================================================

#--------------------------------------------------
if (platform == 'linux-mpi'):
#--------------------------------------------------
   parallel_run = "mpirun -np 4"
   serial_run   = ""
   env = Environment (
      BINPATH     = '#/bin',
      CC          = 'mpicc',	
      CPPDEFINES = ['NO_FREETYPE'],
      CPPFLAGS    = '-Wall -g  -m128bit-long-double',
      CPPPATH     = '#/include',
      CXX         = 'mpiCC',	
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
   env = Environment (
      BINPATH     = '#/bin',
      CC          = 'mpicc',	
      CPPDEFINES = ['NO_FREETYPE'],
      CPPFLAGS    = '-Wall -g  -m128bit-long-double',
      CPPPATH     = '#/include',
      CXX         = 'mpiCC',	
      FORTRAN     = 'gfortran',
      FORTRANLIBS = 'gfortran',
      FORTRANPATH = '#/include',
      LIBPATH     = '#/lib',
   )
#--------------------------------------------------
elif (platform == 'linux-ampi'):
#--------------------------------------------------
   parallel_run = "/home/bordner/charm/charm-6.2.0/bin/charmrun ++p 4 "
   serial_run   = ""
  
   env = Environment(
      BINPATH     = '#/bin',
      CC          = '/home/bordner/charm/charm-6.2.0/bin/charmc -language ampi',
      CPPDEFINES = ['NO_FREETYPE'],
      CPPFLAGS    = '-Wall -g',
      CPPPATH     = '#/include',
      CXX         = '/home/bordner/charm/charm-6.2.0/bin/charmc -language ampi',
      ENV         = os.environ,
      FORTRAN     = 'gfortran',
      FORTRANLIBS = 'gfortran',
      FORTRANPATH = '#/include',
      LIBPATH     = '#/lib',
   )
#--------------------------------------------------
elif (platform == 'triton'):
#--------------------------------------------------
   parallel_run = "/opt/openmpi_pgimx/bin/mpirun -np 4 "
   serial_run   = ""
   env = Environment (
      BINPATH = '#/bin',
      CC      = 'mpicc',	
      CPPDEFINES = ['NO_FREETYPE'],
      CPPFLAGS = '-g -DH5_USE_16_API',
      CPPPATH = ['#/include', '/opt/pgi/hdf5_pgi/include'],
      CXX     = 'mpicxx',	
      ENV = {'PATH' : os.environ['PATH'], 
	'LM_LICENSE_FILE' : os.environ['LM_LICENSE_FILE']},
      FORTRAN = 'mpif90',
      FORTRANPATH = '#/include',
      LIBPATH = ['#/lib', '/opt/pgi/hdf5_pgi/lib'],
      LINKFLAGS = '-pgf90libs',
   )
#--------------------------------------------------
elif (platform == 'ncsa-bd'):
#--------------------------------------------------
   parallel_run = "/opt/openmpi_pgimx/bin/mpirun -np 4 "
   serial_run   = ""
   env = Environment (
      ARFLAGS  = 'r',
      BINPATH = '#/bin',
      CCFLAGS = '-O3 -qhot -q64 -D H5_USE_16_API',
      CC      = 'mpcc',	
      CPPDEFINES = ['NO_FREETYPE'],
#      CPPDEFPREFIX = '-WF,-D',
      CPPPATH = ['/home/bordner/include', '#/include'],
      CXX     = 'mpCC',	
      DEFINES = '',
      FORTRANFLAGS = '-O3 -qhot -q64 -qextname',
      FORTRANLIBS = 'xlf90',
      FORTRAN = 'mpfort',
      LIBPATH = ['#/lib','/home/bordner/lib','/opt/ibmcmp/xlf/13.1/lib64'],
      LINKFLAGS  = '-q64'
   )
elif (platform == 'unknown'):
   print
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
   print "   linux-mpi    Linux with MPI"
   print "   linux-ampi   Linux with AMPI (Charm++)"
   print "   triton:      SCSD Triton Resource"
   print "   ncsa-bd      NCSA Blue Drop"
   print
   print "**********************************************************************"
   print
   sys.exit()
   env = Environment (
      CPPPATH = ['#/include'],
      LIBPATH = ['#/lib'],
      BINPATH = '#/bin'
   )


Export('env')
Export('platform')
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



