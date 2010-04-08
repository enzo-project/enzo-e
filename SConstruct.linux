
TOPDIR  = GetLaunchDir()

print TOPDIR + '/include'

env = Environment (
	CXX     = 'mpiCC',	
	CPPPATH = '#/include',
	LIBPATH = '#/lib',
	BINPATH = '#/bin',
	CPPFLAGS = '-Wall -g'
         )

Export('env')

SConscript('src/SConscript')
SConscript('test/SConscript')



