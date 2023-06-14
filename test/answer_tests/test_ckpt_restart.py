# all tests in this file are intended to test Checkpoint-Restart functionallity
# - at the moment, we are only testing the new-style, scalable machinery
# - however, it would also be good to test the legacy charm++-based machinery
#
# The tests in this file are different from most other tests in this directory
# because:
#   1. They require 2 dependent calls to Enzo-E (the checkpoint run followed by
#      the restart run)
#   2. They aren't actually answer tests
#
# These tests will currently only work with methods that are deterministic &
# bitwise reproducible. So, it won't work for anything like Gravity methods or
# Methods that include a reduction of floating point values (with variable
# ordering)

import os
import tempfile

import pytest

from test_utils.ckpt_restart_testing import run_ckpt_restart_test
from answer_testing import cached_opts, get_symlink_targets

_input_dir = os.path.abspath(
    os.path.abspath(os.path.join(os.path.dirname(__file__), "../../input"))
)

# This is the test for the new-style, scalable checkpoint-restart test
# -

@pytest.mark.parametrize("nominal_input,uses_grackle",
                         [('checkpoint_ppm.in', False),
                          ('checkpoint_grackle.in', True)])
def test_ckpt_restart(nominal_input, uses_grackle):
    nominal_input = os.path.abspath(f'{_input_dir}/Checkpoint/{nominal_input}')

    symlink_srcs = list(get_symlink_targets(grackle_files = uses_grackle))

    # the following context manager class creates a temporary directory,
    # stores the path to it in working_dir, and then always cleans up the
    # contents (even if the program encounters an error) when we exit the
    # context
    #
    # TODO: it would be REALLY nice to have a unified approach (shared across
    #       all tests) for doing this, that can be configured to:
    #       - name the temporary-directory in a deterministic manner after the
    #         current test (maybe let the user specify an empty directory that
    #         all tests get run inside of)
    #       - conditionally decide whether to delete the temporary directory.
    #         One could imagine 3 options:
    #            1. Always delete the temporary directory (current behavior)
    #            2. Only delete the temporary directory on test failure
    #            3. Never delete the temporary directory
    with tempfile.TemporaryDirectory() as working_dir:

        run_ckpt_restart_test(os.path.join(_input_dir, nominal_input),
                              working_dir = working_dir,
                              enzoe_driver = cached_opts().enzoe_driver,
                              nproc = 1, ckpt_cycle = 2, stop_cycle = 4,
                              symlink_srcs = symlink_srcs,
                              sim_name_prefix = None,
                              buffer_outputs_on_disk = False)
