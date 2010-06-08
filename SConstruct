import os
import sys

TOPDIR  = GetLaunchDir()

# Initialize platform

platform='unknown'

platform = ARGUMENTS.get('platform','unknown')

if (platform == 'unknown' and "CELLO_PLATFORM" in os.environ):
     platform = os.environ["CELLO_PLATFORM"]

# Get the SVN code revision numbers

revision_new = os.popen("svn info svn+ssh://client65-88.sdsc.edu/usr/local/svn/cello/trunk | awk '/Revision:/ {print $2}'").read()

revision_current = os.popen("svn info | awk '/Revision:/ {print $2}'").read()

revision_changes = os.popen("svn status | grep -v '?' | wc -l").read()

if (revision_changes != "0"):
   print "\nWARNING: Working directory has local modifications!\n"

if (revision_new != revision_current):
   print "WARNING: Working directory is not up-to-date with the repository!"


# Initialize environment according to platform

#--------------------------------------------------
if (platform == 'linux-mpi'):
#--------------------------------------------------
   parallel_run = "mpirun -np 4"
   serial_run   = ""
   env = Environment (
      CXX         = 'mpiCC',	
      CPPFLAGS    = '-Wall -g  -m128bit-long-double',
#      CPPFLAGS    = '-Wall -g  -m96bit-long-double',
      CPPPATH     = '#/include',
      FORTRANPATH = '#/include',
      FORTRAN     = 'gfortran',
      FORTRANLIBS = 'gfortran',
      LIBPATH     = '#/lib',
      BINPATH     = '#/bin',
   )
#--------------------------------------------------
elif (platform == 'linux-mpi-valgrind'):
#--------------------------------------------------
   parallel_run = "mpirun -np 4 valgrind"
   serial_run   = "valgrind "
   env = Environment (
      CXX         = 'mpiCC',	
      CPPFLAGS    = '-Wall -g  -m128bit-long-double',
#      CPPFLAGS    = '-Wall -g  -m96bit-long-double',
      CPPPATH     = '#/include',
      FORTRANPATH = '#/include',
      FORTRAN     = 'gfortran',
      FORTRANLIBS = 'gfortran',
      LIBPATH     = '#/lib',
      BINPATH     = '#/bin',
   )
#--------------------------------------------------
elif (platform == 'linux-ampi'):
#--------------------------------------------------
   parallel_run = "/home/bordner/charm/charm-6.2.0/bin/charmrun ++p 4 "
   serial_run   = ""
  
   env = Environment(
      ENV         = os.environ,
      CXX         = '/home/bordner/charm/charm-6.2.0/bin/charmc -language ampi',
      CPPFLAGS    = '-Wall -g',
      CPPPATH     = '#/include',
      FORTRANPATH = '#/include',
      FORTRAN     = 'gfortran',
      FORTRANLIBS = 'gfortran',
      LIBPATH     = '#/lib',
      BINPATH     = '#/bin',
   )
#--------------------------------------------------
elif (platform == 'triton'):
#--------------------------------------------------
   parallel_run = "/opt/openmpi_pgimx/bin/mpirun -np 4 "
   serial_run   = ""
   env = Environment (
      ENV = {'PATH' : os.environ['PATH'],
	       'LM_LICENSE_FILE' : os.environ['LM_LICENSE_FILE']},
      CXX     = 'mpicxx',	
      CPPPATH = ['#/include', '/opt/pgi/hdf5_pgi/include'],
      LIBPATH = ['#/lib', '/opt/pgi/hdf5_pgi/lib'],
      BINPATH = '#/bin',
      CPPFLAGS = '-g -DH5_USE_16_API',
      FORTRAN = 'mpif90',
      FORTRANPATH = '#/include',
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



