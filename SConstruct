
TOPDIR  = GetLaunchDir()

print TOPDIR + '/include'

env = Environment (
	CXX     = 'mpiCC',	
	CPPPATH = '#/include',
	LIBPATH = '#/lib',
	BINPATH = '#/bin'
         )

Export('env')

SConscript('bin/SConscript')
SConscript('src/SConscript')



