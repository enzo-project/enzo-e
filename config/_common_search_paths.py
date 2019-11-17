import os

def charm_path_search(home):
	# function that searches for charm_path in locations common across
	# multiple architectures
	charm_path = os.getenv('CHARM_HOME')

	if (charm_path is None) and (home is not None):
		if os.path.isdir(home + '/Charm/charm'):
			charm_path = home + '/Charm/charm'
		elif os.path.isdir(home + '/charm'):
			charm_path = home + '/charm'
		elif os.path.isdir(home + '/local/charm'):
			charm_path = home + '/local/charm'
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
			raise Exception('Charm++ was not found.	 Try setting the CHARM_HOME environment variable.')
	return charm_path

def grackle_path_search(home):
	# function that searches for grackle_path in locations common across
	# multiple architectures
	grackle_path = os.getenv('GRACKLE_HOME')
	if (grackle_path is None) and (home is not None):
		if os.path.isdir(home + '/Grackle/src/clib'):
			grackle_path = home + '/Grackle/src/clib'
		elif os.path.isdir(home + '/grackle/src/clib'):
			grackle_path = home + '/grackle/src/clib'
		elif os.path.isdir(home + '/local/grackle/src/clib'):
			grackle_path = home + '/local/grackle/src/clib'
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
	return grackle_path
