#!/usr/bin/python3

import argparse
import os
import os.path
import shutil
import sys

_test_utils_path = os.path.abspath(os.path.join(
    os.path.dirname(__file__), "../test/answer_tests"
))
sys.path.append(_test_utils_path)

from test_utils.enzoe_driver import EnzoEDriver
from test_utils.ckpt_restart_testing import run_ckpt_restart_test

def main(args):
    
    _charmpath = None
    if args.charm is not None:
        _charmpath = os.path.abspath(args.charm)

    enzoe_driver = EnzoEDriver(enzoe_path = os.path.abspath(args.enzoe),
                               charmrun_path = _charmpath)

    # determine the directory where we will run the tests & cleanup policy
    if args.test_dir is None:
        from tempfile import mkdtemp
        work_dir = mkdtemp(prefix = 'ckpt-restart-scratch-', dir = './')
        cleanup_policy = 'always' if (args.cleanup is None) else args.cleanup
    else:
        work_dir = args.test_dir
        cleanup_policy = 'never' if (args.cleanup is None) else args.cleanup

    # check that work_dir exists & is empty
    try:
        contents = [os.path.join(work_dir,e) for e in os.listdir(work_dir)]
    except FileNotFoundError:
        raise RuntimeError(f"TEST_DIR, '{work_dir}', doesn't exist") from None
    for path in contents: # docs warn against using os.scandir with deletion
        if not args.clobber:
            raise RuntimeError(f"TEST_DIR, '{work_dir}', must be empty")
        elif os.path.islink(path) or os.path.isfile(path):
            os.remove(path)
        else:
            shutil.rmtree(path)

    # now, run the test:
    test_complete = False
    try:
        run_ckpt_restart_test(nominal_input = args.input,
                              working_dir = work_dir,
                              enzoe_driver = enzoe_driver, nproc = 1,
                              ckpt_cycle = args.restart_cycle,
                              stop_cycle = args.stop_cycle,
                              symlink_srcs = args.symlinks,
                              sim_name_prefix = None)
        test_complete = True
    finally:
        perform_cleanup = ( (cleanup_policy == 'always') or
                            (test_complete and (cleanup_policy == 'success')) )

        print(f"Test {['failed', 'passed'][test_complete]}: " +
              ["NOT cleaning up", "cleaning up"][perform_cleanup])
        if perform_cleanup:
            shutil.rmtree(work_dir)


# Setup the parser for the CLI
_description = """\
Perform a checkpoint-restart test.

The test performs the following procedure:
  1. Parameters are taken from the specified Cello/Enzo-E input file and 
     modified (without affecting the original) copied such that a simulation
     writes a checkpoint output to disk at RESTART_CYCLE, writes a data output
     to disk including all fields at STOP_CYCLE-1 and runs until STOP_CYCLE.
  2. The simulation is run using the modified parameters from within
     TEST_DIR/ckpt_run.
  3. The simulation is then restarted while running from within
     TEST_DIR/restart_run.
  4. Finally, the program ensures that the output files written at STOP_CYCLE
     contain identical data before and after the restart.

Note: This program executes Cello/Enzo-E in a temporary directory, so that
multiple simultaneous executions of this test should not interfere with each
other.
"""

parser = argparse.ArgumentParser(
    description = _description,
    formatter_class=argparse.RawDescriptionHelpFormatter
)
parser.add_argument('-i','--input', required = True,
                    help = ('Required option used to specify the template '
                            'input file. This can\'t have any scalar/logical '
                            'expressions or append lists to parameters.'))
parser.add_argument('--test-dir', default = None,
                    help = ('Name of an existing directory that will contain '
                            'all scratch-outputs produced during the test. If '
                            'not specified, a temporary directory is made'))
parser.add_argument('-r','--restart_cycle', type = int, default = 2,
                    help = ('The cycle from which the simulation restarts. '
                            'Default is 2.'))
parser.add_argument('-s','--stop_cycle', type = int, default = 4,
                    help = ('The cycle at which data outputs are compared. '
                            'Default is 4.'))
parser.add_argument('--enzoe', required = True, help = 'path to enzo-e binary')
parser.add_argument('--charm', default = None,
                    help = ('the path to charmrun. If not specified, the '
                            'simulation will be run serially.'))
parser.add_argument('--symlinks', nargs = '+', default = [],
                    help = ('paths to files/directories for which we will '
                            'create symlinks in each run directory'))

parser.add_argument('--clobber', action = 'store_true', dest = 'clobber',
                    help = ('Remove all contents from TEST_DIR (if it is user '
                            'specified), before doing anything else. By '
                            'default, the program fails if TEST_DIR isn\'t '
                            'empty.'))
parser.add_argument('--cleanup', choices = ['always', 'never', 'success'],
                    default = None,
                    help = ('Specifies when to delete TEST_DIR. When the user '
                            'specifies TEST_DIR, the default is "never". '
                            'Otherwise default is "always".'))

if __name__ == '__main__':
    main(parser.parse_args())
