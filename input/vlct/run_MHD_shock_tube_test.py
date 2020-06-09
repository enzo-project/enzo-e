#!/bin/python

# runs VLCT MHD rj2a Shock Tube Tests
# - This should probably be transitioned to testing framework once it's created
# - This script expects to be called from the root level of the repository
#   OR at the same level where its defined
#
# These tests draw loose inspiration from shock tube tests used by Athena++.
# Specifically this script:
#   1.) Checks that the L1-norm of the error for the from Ryu & Jones (1995)
#       shock tube problem when the problem evolves along all 3 axes.
#
# We should probably try small convergence tests like Athena++

import os.path
import sys
import shutil
from functools import partial

import numpy as np

from testing_utils import \
    prep_cur_dir, standard_analyze, EnzoEWrapper, CalcTableL1Norm

def run_tests(executable):

    temp = 'input/vlct/MHD_shock_tube/method_vlct_{:s}_rj2a_N{:d}.in'
    call_test = EnzoEWrapper(executable,temp)

    # can't actually run the serial vs parallel test from this script (not sure
    # how to get around request for password)
    call_test("x",256)
    call_test("y",256)
    call_test("z",256)

def analyze_shock(ref_val, axis, target_template, name_template,
                  l1_functor, reverse = False, std_dev = False,
                  offset = 0., bkg_vel=(), verbose = False):
    """
    Helper function for simplifying permutable axis-aligned shock tube test
    analysis.

    Notes
    -----
    The offset and bkg_vel kwargs have been introduced for test problems in
    which shock tubes have a background velocity in the frame where the mesh is
    motionless (used for testing dual energy).
    """
    axis = axis.lower()
    assert axis in 'xyz'
    assert len(bkg_vel) in [0,3]
    functor_kwargs = dict(axis = axis, std_dev = std_dev, reverse = reverse,
                          bkg_velocity = bkg_vel, permute = (axis != 'x'))
    if offset != 0.:
        functor_kwargs['offset_soln'] = offset

    return standard_analyze(ref_val, l1_functor, target_template.format(axis),
                            name_template.format(axis),
                            functor_kwargs = functor_kwargs,
                            exact = std_dev, verbose = verbose)

def analyze_tests():
    # define the functor for evaluating the norm of the L1 error vector
    ref_table = "input/vlct/MHD_shock_tube/rj2a_shock_tube_t0.2_res256.csv"
    l1_func= CalcTableL1Norm("tools/l1_error_norm.py",
                             ["density","velocity_x","velocity_y","velocity_z",
                              "pressure","bfield_x","bfield_y","bfield_z"],
                             default_ref_table = ref_table)

    err_compare = partial(analyze_shock,
                          target_template = "method_vlct-1-{:s}_rj2a_N256_0.2",
                          name_template = "{}-axis rj2a shock tube N=256",
                          l1_functor = l1_func)

    r = []
    # check the average L1-Norm along the active axis (averaged over multiple
    # slices along the active axis)
    # check that the standard deviation of the L1-Norms computed along the
    # active axis (There should be no definitely be no differences if we only
    # use 1 block - if we use more than one block, it's unclear to me if it's
    # ok to have round-off errors)

    r.append(err_compare(0.012524502240006011,"x"))
    r.append(err_compare(0.0, "x", std_dev=True))
    r.append(err_compare(0.012524502240005972,"y"))
    r.append(err_compare(0.0, "y", std_dev=True))
    r.append(err_compare(0.012524502240005921,"z"))
    r.append(err_compare(0.0, "z", std_dev=True))

    n_passed = np.sum(r)
    n_tests = len(r)
    print("{:d} Tests passed out of {:d} Tests.".format(n_passed,n_tests))

    return n_passed == n_tests

def cleanup():

    dir_names = ["method_vlct-1-x_rj2a_N256_0.2",
                 "method_vlct-1-y_rj2a_N256_0.2",
                 "method_vlct-1-z_rj2a_N256_0.2"]
    for dir_name in dir_names:
        if os.path.isdir(dir_name):
            shutil.rmtree(dir_name)

if __name__ == '__main__':

    executable = 'bin/enzo-p'

    # this script can either be called from the base repository or from
    # the subdirectory: input/vlct
    prep_cur_dir(executable)

    # run the tests
    run_tests(executable)

    # analyze the tests
    tests_passed = analyze_tests()

    # cleanup the tests
    cleanup()

    if tests_passed:
        sys.exit(0)
    else:
        sys.exit(3)
