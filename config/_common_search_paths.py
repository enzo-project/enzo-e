import os

def _check_dirs(search_paths):
        for search_path in search_paths:
                if os.path.isdir(search_path):
                        return path
        return None

def charm_path_search(home):
	# function that searches for charm_path in locations common across
	# multiple architectures
	charm_path = os.getenv('CHARM_HOME')

	if (charm_path is None) and (home is not None):
                temp = ['Charm/charm', 'charm', 'local/charm', 'src/charm',
                        'source/charm']
                paths = [os.path.join(home, elem) for elem in temp]
                charm_path = _check_dirs(paths)
	if charm_path is None:
                charm_path = _check_dirs(['/usr/local/charm', '/opt/charm'])
                if charm_path is None:
			raise Exception('Charm++ was not found.	Try setting '
                                        'the CHARM_HOME environment variable.')
	return charm_path

def grackle_path_search(home):
	# function that searches for grackle_path in locations common across
	# multiple architectures
	grackle_path = os.getenv('GRACKLE_HOME')
	if (grackle_path is None) and (home is not None):
                temp = ['Grackle', 'grackle', 'local/grackle',
                        'src/grackle', 'source/grackle', 'Software/Grackle']
                paths = [os.path.join(home, elem, 'src/clib') for elem in temp]
                grackle_path = _check_dirs(paths)
	if grackle_path is None:
                grackle_path = _check_dirs(['/usr/local/grackle/src/clib',
                                            '/opt/grackle/src/clib'])
	return grackle_path
