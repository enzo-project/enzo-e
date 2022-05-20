#!/bin/python

# runs VLCT HD inclined stable Jeans Wave Tests
# - This should probably be transitioned to testing framework once it's created
# - This script expects to be called from the root level of the repository
#   OR at the same level where its defined
#
# These tests draw loose inspiration from Linear Wave tests used by Athena++.
# Specifically, this script:
#   1.) Checks if that L1-norm of the error is appropriate for the inclined
#       stable jeans wave for N=16 and N=32

import os
import os.path
import shutil
import sys

from functools import partial

import numpy as np

from testing_utils import \
    CalcSimL1Norm, EnzoEWrapper, isclose, standard_analyze, testing_context

from run_MHD_linear_wave_test import \
    analyze_linwave, identical_l1_error_linwave

def run_tests(executable):
    temp = 'input/vlct/jeans_wave/vlct_stable_jeansN{:d}.in'
    call_test = EnzoEWrapper(executable,temp)

    call_test(16)
    call_test(32)

def analyze_tests():
    # define the functor for evaluating the norm of the L1 error vector
    l1_func = CalcSimL1Norm(["density","velocity_x","velocity_y","velocity_z",
                             "pressure",])
    template = "method_vlct-cg-1-inclined-jeansN{res}_{time:.4f}"

    verbose = False
    err_compare = partial(analyze_linwave, target_template = template,
                          name_template = 'stable-jeans wave N={res}',
                          l1_functor = l1_func,
                          final_times = {'stable-jeans' : 1.0},
                          verbose = verbose)
    r = []
    # first let's run the l1-norm of each value for 2 sizes to make sure we get
    # convergence (originally we were looking at the slope of convergence but
    # that allows the errors to be larger)

    r.append(err_compare(1.2316886043695045e-07, 'stable-jeans', 1,16))
    r.append(err_compare(1.903655646665732e-08, 'stable-jeans', 1,32))

    n_passed = np.sum(r)
    n_tests = len(r)
    success = (n_passed == n_tests)
    print("{:d} Tests passed out of {:d} Tests.".format(n_passed,n_tests))
    return success

def cleanup():

    dir_names = ["method_vlct-cg-1-inclined-jeansN16_0.0000",
                 "method_vlct-cg-1-inclined-jeansN16_1.0000",
                 "method_vlct-cg-1-inclined-jeansN32_0.0000",
                 "method_vlct-cg-1-inclined-jeansN32_1.0000"]

    for dir_name in dir_names:
        if os.path.isdir(dir_name):
            shutil.rmtree(dir_name)

if __name__ == '__main__':

    executable = os.environ.get('ENZOE_BIN', 'bin/enzo-e')

    with testing_context():
        # run the tests
        tests_complete = run_tests(executable)

        # analyze the tests
        tests_passed = analyze_tests()

        # cleanup the tests
        cleanup()

    if tests_passed:
        sys.exit(0)
    else:
        sys.exit(3)
