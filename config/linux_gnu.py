import os

is_arch_valid = 1

flags_arch = '-Wall -O3 -g'
#flags_arch = '-Wall -pg -fprofile-arcs -ftest-coverage'
#flags_arch = '-Wall -pg'
#flags_arch = '-Wall -g'
flags_link_charm = ' -rdynamic' # required for backtraces

cc  = 'gcc '
f90 = 'gfortran'

flags_prec_single = '-fdefault-real-4 -fdefault-double-8'
flags_prec_double = '-fdefault-real-8 -fdefault-double-8'

libpath_fortran = '.'
libs_fortran    = ['gfortran']

home = os.getenv('HOME')

charm_path = os.getenv('CHARM_HOME')
if charm_path is None:
	if home is not None:
		if os.path.isdir(home + '/Charm/charm'):
			charm_path = home + '/Charm/charm'
		elif os.path.isdir(home + '/charm'):
			charm_path = home + '/charm'
		elif os.path.isdir(home + '/src/charm'):
			charm_path = home + '/src/charm'
		elif os.path.isdir(home + '/source/charm'):
			charm_path = home + '/source/charm'
	if charm_path is None:
		if os.path.isdir('/usr/local/charm'):
			charm_path = '/usr/local/charm'
		elif os.path.isdir('/opt/charm'):
			charm_path = '/opt/charm'
		else:
			raise Exception('Charm++ was not found.  Try setting the CHARM_HOME environment variable.')

papi_path    = '/usr/local'
hdf5_path    = '/usr'
png_path     = '/lib/x86_64-linux-gnu'

grackle_path = os.getenv('GRACKLE_HOME')
if grackle_path is None:
	if home is not None:
		if os.path.isdir(home + '/Grackle/src/clib'):
			grackle_path = home + '/Grackle/src/clib'
		elif os.path.isdir(home + '/grackle/src/clib'):
			grackle_path = home + '/grackle/src/clib'
		elif os.path.isdir(home + '/src/grackle/src/clib'):
			grackle_path = home + '/src/grackle/src/clib'
		elif os.path.isdir(home + '/source/grackle/src/clib'):
			grackle_path = home + '/source/grackle/src/clib'
		elif os.path.isdir(home + '/Software/Grackle/src/clib'):
			grackle_path = home + '/Software/Grackle/src/clib'
	if grackle_path is None:
		if os.path.isdir('/usr/local/grackle/src/clib'):
			grackle_path = '/usr/local/grackle/src/clib'
		elif os.path.isdir('/opt/grackle/src/clib'):
			grackle_path = '/opt/grackle/src/clib'
