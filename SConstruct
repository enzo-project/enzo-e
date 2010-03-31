
TOPDIR  = GetLaunchDir()

print TOPDIR + '/include'

env = Environment (
	CXX     = 'mpiCC',	
	CPPPATH = '#/include',
	LIBPATH = '#/lib',
	BINPATH = '#/bin'
         )

Export('env')

SConscript('src/SConscript')

