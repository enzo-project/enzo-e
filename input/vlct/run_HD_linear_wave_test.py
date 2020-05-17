#!/bin/python

# runs VLCT HD Linear Wave Tests
# - This should probably be transitioned to testing framework once it's created
# - This script expects to be called from the root level of the repository
#   OR at the same level where its defined
#
# These tests draw loose inspiration from Linear Wave tests used by Athena++.
# Specifically, this script:
#   1.) Checks if that L1-norm of the error is appropriate for the inclined
#       fast wave, alfven wave, slow wave, and entropy for N=16 and N=32
#       [Unlike Athena++, we are compare L1-norms to precise expected values - 
#        therefore it probably isn't necessary to use 2 different resolutions]
#   2.) Check that the L1-norm error is the same for left propagating and right
#       propagating fast waves

import os.path
import shutil
import sys

from functools import partial

import numpy as np

from testing_utils import \
  CalcSimL1Norm, EnzoEWrapper, isclose, standard_analyze, prep_cur_dir

from run_linear_wave_test import analyze_linwave, identical_l1_error_linwave

def run_tests(executable):
    temp = 'input/vlct/HD_linear_wave/method_vlct_{}N{:d}{}.in'
    call_test = EnzoEWrapper(executable,temp)

    call_test("sound",16,'')
    call_test("sound",32,'')
    call_test("sound",32,'_left')

    call_test("hd_entropy",16,'')
    call_test("hd_entropy",32,'')

    call_test("hd_transv_entropy_v1",16,'')
    call_test("hd_transv_entropy_v1",32,'')

    call_test("hd_transv_entropy_v2",16,'')
    call_test("hd_transv_entropy_v2",32,'')

def analyze_tests():
    # define the functor for evaluating the norm of the L1 error vector
    l1_func = CalcSimL1Norm("tools/l1_error_norm.py",
                            ["density","velocity_x","velocity_y","velocity_z",
                             "pressure","bfield_x","bfield_y","bfield_z"])

    # define the template for the directory holding the simulation data
    template = "method_vlct-{nblocks}-{wave_name}N{res}{direction}_{time:.1f}"

    final_times = {"sound" : 1.0, "hd_entropy" : 1.0,
                   "hd_transv_entropy_v1" : 1.0, "hd_transv_entropy_v2" : 1.0}

    verbose = False
    err_compare = partial(analyze_linwave, target_template = template,
                          name_template = '{wave_name} wave N={res}',
                          l1_functor = l1_func, final_times = final_times,
                          verbose = verbose)

    r = []
    # first let's run the l1-norm of each value for 2 sizes to make sure we get
    # convergence (originally we were looking at the slope of convergence but
    # that allows the errors to be larger)

    r.append(err_compare(1.3704437788774914e-07, 'sound',                1,16))
    r.append(err_compare(2.348494680480983e-08,  'sound',                1,32))

    r.append(err_compare(8.736217559091042e-08,  'hd_entropy',           1,16))
    r.append(err_compare(2.2383944940640457e-08, 'hd_entropy',           1,32))

    r.append(err_compare(7.671004076383043e-08,  'hd_transv_entropy_v1', 1,16))
    r.append(err_compare(1.6725791096099574e-08, 'hd_transv_entropy_v1', 1,32))

    r.append(err_compare(8.64073017675657e-08,   'hd_transv_entropy_v2', 1,16))
    r.append(err_compare(1.6960029746329912e-08, 'hd_transv_entropy_v2', 1,32))

    # Check error between left and right propagating waves. As in Athena++ only
    # check to 6 decimal places
    r.append(identical_l1_error_linwave(1, 'sound', 32, l1_func, template, 1.0,
                                        prec=6, verbose = verbose))
    n_passed = np.sum(r)
    n_tests = len(r)
    success = (n_passed == n_tests)
    print("{:d} Tests passed out of {:d} Tests.".format(n_passed,n_tests))
    return success

def cleanup():

    dir_names = ["method_vlct-1-soundN16_0.0",
                 "method_vlct-1-soundN16_1.0",
                 "method_vlct-1-soundN32_0.0",
                 "method_vlct-1-soundN32_1.0",
                 "method_vlct-1-soundN32-left_0.0",
                 "method_vlct-1-soundN32-left_1.0",
                 "method_vlct-1-hd_entropyN16_0.0",
                 "method_vlct-1-hd_entropyN16_1.0",
                 "method_vlct-1-hd_entropyN32_0.0",
                 "method_vlct-1-hd_entropyN32_1.0",
                 "method_vlct-1-hd_transv_entropy_v1N16_0.0",
                 "method_vlct-1-hd_transv_entropy_v1N16_1.0",
                 "method_vlct-1-hd_transv_entropy_v1N32_0.0",
                 "method_vlct-1-hd_transv_entropy_v1N32_1.0",
                 "method_vlct-1-hd_transv_entropy_v2N16_0.0",
                 "method_vlct-1-hd_transv_entropy_v2N16_1.0",
                 "method_vlct-1-hd_transv_entropy_v2N32_0.0",
                 "method_vlct-1-hd_transv_entropy_v2N32_1.0"]

    for dir_name in dir_names:
        if os.path.isdir(dir_name):
            shutil.rmtree(dir_name)

if __name__ == '__main__':

    executable = 'bin/enzo-p'

    # this script can either be called from the base repository or from
    # the subdirectory: input/vlct
    prep_cur_dir(executable)

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
