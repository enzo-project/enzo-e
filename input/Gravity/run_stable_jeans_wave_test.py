#!/bin/python

# runs VLCT HD and PPM HD inclined stable Jeans Wave Tests
# - This should probably be transitioned to testing framework once it's created
# - This script expects to be called from the root level of the repository
#   OR at the same level where its defined
#
# These tests draw loose inspiration from Linear Wave tests used by Athena++.
# Specifically, this script: Checks if that L1-norm of the error is appropriate
# for the inclined stable jeans wave for N=16 and N=32
#
# In the future, we may want to split the PPM and VL+CT tests

import argparse
import os.path
import shutil
import sys

from functools import partial

import numpy as np

# import testing utilities defined for VL+CT tests (this approach is very hacky
# - we really need to revisit this in the future!)
_LOCAL_DIR = os.path.dirname(os.path.realpath(__file__))
_VLCT_DIR = os.path.join(_LOCAL_DIR, "../vlct")
if os.path.isdir(_VLCT_DIR):
    sys.path.insert(0, _VLCT_DIR)
    from testing_utils import \
        CalcSimL1Norm, EnzoEWrapper, isclose, standard_analyze, testing_context
    from run_MHD_linear_wave_test import \
        analyze_linwave, identical_l1_error_linwave
else:
    raise RuntimeError(f"expected VL+CT tests to be defined in {_VLCT_DIR}, "
                       "but that that directory does not exist")

def run_tests(executable):

    for solver in ['ppm', 'vlct']:

        temp = 'input/Gravity/jeans_wave/' + solver + '_stable_jeansN{:d}.in'
        call_test = EnzoEWrapper(executable,temp)

        call_test(16)
        call_test(32)

def analyze_tests():
    # define the functor for evaluating the norm of the L1 error vector
    l1_func = CalcSimL1Norm(["density","velocity_x","velocity_y","velocity_z",
                             "pressure",])

    verbose = False
    err_cmp_vlct = partial(
        analyze_linwave,
        target_template = "method_vlct-cg-1-inclined-jeansN{res}_{time:.4f}",
        name_template = 'vlct-stable-jeans wave N={res}', l1_functor = l1_func,
        final_times = {'vlct-stable-jeans' : 1.1547}, verbose = verbose
    )
    err_cmp_ppm = partial(
        analyze_linwave,
        target_template = "method_ppm-cg-1-inclined-jeansN{res}_{time:.4f}",
        name_template = 'ppm-stable-jeans wave N={res}', l1_functor = l1_func,
        final_times = {'ppm-stable-jeans' : 1.1547}, verbose = verbose
    )

    r = []
    # first let's run the l1-norm of each value for 2 sizes to make sure we get
    # convergence (originally we were looking at the slope of convergence but
    # that allows the errors to be larger)

    r.append(err_cmp_vlct(5.104682632606632e-07, 'vlct-stable-jeans', 1,16))
    r.append(err_cmp_vlct(2.2411686465108514e-07, 'vlct-stable-jeans', 1,32))
    r.append(err_cmp_ppm(1.6245498643250925e-07, 'ppm-stable-jeans', 1,16))
    r.append(err_cmp_ppm(7.379045545039628e-08, 'ppm-stable-jeans', 1,32))

    n_passed = np.sum(r)
    n_tests = len(r)
    success = (n_passed == n_tests)
    print("{:d} Tests passed out of {:d} Tests.".format(n_passed,n_tests))
    return success

def cleanup():

    dir_names = ["method_vlct-cg-1-inclined-jeansN16_0.0000",
                 "method_vlct-cg-1-inclined-jeansN16_1.0000",
                 "method_vlct-cg-1-inclined-jeansN32_0.0000",
                 "method_vlct-cg-1-inclined-jeansN32_1.0000",
                 "method_ppm-cg-1-inclined-jeansN16_0.0000",
                 "method_ppm-cg-1-inclined-jeansN16_1.0000",
                 "method_ppm-cg-1-inclined-jeansN32_0.0000",
                 "method_ppm-cg-1-inclined-jeansN32_1.0000"]

    for dir_name in dir_names:
        if os.path.isdir(dir_name):
            shutil.rmtree(dir_name)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--launch_cmd', required=True,type=str)
    args = parser.parse_args()

    with testing_context():
        # run the tests
        tests_complete = run_tests(args.launch_cmd)

        # analyze the tests
        tests_passed = analyze_tests()

        # cleanup the tests
        cleanup()

    if tests_passed:
        sys.exit(0)
    else:
        sys.exit(3)
