#!/bin/python

# runs VLCT MHD rj2a Shock Tube Tests
# - This should probably be transitioned to testing framework once it's created
# - This script expects to be called from the root level of the repository
#   OR at the same level where its defined
#
# These tests draw loose inspiration from Linear Wave tests used by Athena++.
# Specifically this script:
#   1.) Checks that the L1-norm of the error for the from Ryu & Jones (1995)
#       shock tube problem when the problem evolves along all 3 axes.
#
# After the bug get's fixed with using periodic boundary conditions alongside
# outflow conditions, we should probably try small convergence tests like
# Athena++

import os
import os.path
import numpy as np
import shutil
import subprocess
import math

from run_linear_wave_test import isclose, prep_cur_dir

# this executes things in standalone mode
_executable = 'bin/enzo-p'
l1_norm_table_template = "python tools/l1_error_norm.py table{:s} {:s} {:s}"

data_dir_template = "method_vlct-1-{:s}_rj2a_N{:d}_{:.1f}"

def call_test(axis,res):
    input_file_temp = 'input/vlct/shock_tube/method_vlct_{:s}_rj2a_N{:d}.in'
    input_file = input_file_temp.format(axis, res)

    command = _executable + ' ' + input_file
    subprocess.call(command,shell=True)

def run_tests():

    # can't actually run the serial vs parallel test from this script (not sure
    # how to get around request for password)
    call_test("x",256)
    call_test("y",256)
    call_test("z",256)

def calc_table_l1_norm(axis,table_path,target_path, std_dev = False,
                       permute = False):
    command = l1_norm_table_template.format(axis, table_path, target_path)
    if std_dev:
        command = command + " --std"
    if permute:
        command = command + " --permute"

    return float(subprocess.check_output(command,shell=True))

def standard_rj2a_l1_analyze(axis,res,ref_l1_norm, std_dev = False,
                             exact = False):
    ref_table = "input/vlct/shock_tube/rj2a_shock_tube_t0.2_res256.csv"
    target_path = data_dir_template.format(axis,res,0.2)
    norm = calc_table_l1_norm(axis,ref_table,target_path,std_dev=std_dev,
                              permute = (axis != 'x'))

    if exact and std_dev and norm != ref_l1_norm:
        message_temp = ("Standard Deviation of norm of L1 error of {}-axis "
                        "rj2a shock tube N={:d} isn't correct\n{:s} {:s}")
        print(message_temp.format(axis, res, repr(norm), repr(ref_l1_norm)))
        return False
    if exact != std_dev:
        raise NotImplementedError()
    if (not exact) and (not isclose(norm, ref_l1_norm, abs_tol = True)):
        message_temp = ("L1 error of {}-axis rj2a shock tube N={:d} isn't "
                        "correct\n{:s} {:s}")
        print(message_temp.format(axis, res, repr(norm), repr(ref_l1_norm)))
        return False
    return True

def analyze_tests():

    # the values are taken from a local run

    r = []
    # check the average L1-Norm along the active axis (averaged over multiple
    # slices along the active axis)
    # check that the standard deviation of the L1-Norms computed along the
    # active axis (There should be no definitely be no differences if we only
    # use 1 block - if we use more than one block, it's unclear to me if it's
    # ok to have round-off errors)

    r.append(standard_rj2a_l1_analyze("x", 256, 0.012524558844892638))
    r.append(standard_rj2a_l1_analyze("x", 256, 0.0, std_dev=True, exact=True))
    r.append(standard_rj2a_l1_analyze("y", 256, 0.012524558844892614))
    r.append(standard_rj2a_l1_analyze("y", 256, 0.0, std_dev=True, exact=True))
    r.append(standard_rj2a_l1_analyze("z", 256, 0.012524558844892873))
    r.append(standard_rj2a_l1_analyze("z", 256, 0.0, std_dev=True, exact=True))

    n_passed = np.sum(r)
    n_tests = len(r)
    print("{:d} Tests passed out of {:d} Tests.".format(n_passed,n_tests))

def cleanup():

    dir_names = ["method_vlct-1-x_rj2a_N256_0.2",
                 "method_vlct-1-y_rj2a_N256_0.2",
                 "method_vlct-1-z_rj2a_N256_0.2"]
    for dir_name in dir_names:
        if os.path.isdir(dir_name):
            shutil.rmtree(dir_name)

if __name__ == '__main__':

    # this script can either be called from the base repository or from
    # the subdirectory: input/vlct
    prep_cur_dir()

    # run the tests
    run_tests()

    # analyze the tests
    analyze_tests()

    # cleanup the tests
    cleanup()
