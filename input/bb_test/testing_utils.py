# Modified version of input/vlct/testing_utils.py. D

# Defines a context manager used by run_bb_test.py

from contextlib import contextmanager
import os
import os.path

try:
    basestring
except NameError:
    basestring = str

import numpy as np

# determine Enzo-E's root directory
if "/input/bb_test" ==  os.path.dirname(os.path.abspath(__file__))[-14:]:
    # this will work even if this file is imported by modifying sys.path
    _ENZOE_ROOT_DIR = os.path.dirname(os.path.abspath(__file__))[:-14]
else:
    raise RuntimeError("run_bb_test.py has been moved. "
                       "Please update the logic for identifying the Enzo-E "
                       "root directory")

@contextmanager
def testing_context(require_enzoe_inputdir = True):
    """
    Context manager to help prepare the current directory for running tests.

    This mainly checks to see whether `./input` is a valid path
      - if it doesn't exist, this creates a symlink to the input directory of
        enzo-e. Upon exitting this context, the symlink is deleted.
      - if `./input` already exists and `require_enzoe_inputdir` is True, this
        ensures that the `./input` is the input directory in the root directory
        of enzo-e or is a symlink to that directory
    """

    path = 'input'

    cleanup = False
    if os.path.isfile(path):  # path is allowed to be a symlink to a dir
        raise RuntimeError('./' + path + ' is a path to a file.')
    elif os.path.isdir(path): # path is allowed to be a symlink to a dir
        realpath = os.path.abspath(os.path.realpath(path))
        expected = os.path.abspath(os.path.join(_ENZOE_ROOT_DIR, 'input'))
        if require_enzoe_inputdir and (realpath != expected):
            raise RuntimeError('./' + path + " doesn't refer to " + expected)
    elif os.path.islink(path):
        raise RuntimeError('./' + path + ' is a broken link.')
    else: # make a symlink to {_ENZOE_ROOT_DIR}/input
        cleanup = True
        os.symlink(src = os.path.join(_ENZOE_ROOT_DIR, path),
                   dst = path, target_is_directory = True)

    try:
        yield None
    finally:
        if cleanup:
            os.unlink(path)
