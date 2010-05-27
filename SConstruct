import os
import sys

TOPDIR  = GetLaunchDir()

# Initialize platform

platform='unknown'

platform = ARGUMENTS.get('platform','unknown')

if (platform == 'unknown' and "CELLO_PLATFORM" in os.environ):
     platform = os.environ["CELLO_PLATFORM"]


# Initialize environment according to platform

#--------------------------------------------------
if (platform == 'linux-mpi'):
#--------------------------------------------------
   env = Environment (
      CXX         = 'mpiCC',	
      CPPFLAGS    = '-Wall -g',
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
SConscript('src/SConscript')
SConscript('test/SConscript')



