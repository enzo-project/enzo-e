#!/bin/python

# runs VLCT MHD Passive Advection sound wave test
# - This should probably be transitioned to testing framework once it's created
# - This script expects to be called from the root level of the repository
#   OR at the same level where its defined
#
# Specifically, this script:
#   1.) Checks if that L1-norm of the error is appropriate for a sound wave
#       that has advected the a sound wave with a passively advected scalar
#       A sin wave is set up in passive scale (in terms of mass fraction) that
#       pi/2 radians out of phase with the rest of the quantities


import os
import os.path
import shutil
import subprocess
import sys

import numpy as np

from testing_utils import EnzoEWrapper, CalcSimL1Norm, isclose, prep_cur_dir

data_dir_template = "method_vlct-1-{:s}_passive_soundN16-{:.1f}"

calc_l1_norm = CalcSimL1Norm("tools/l1_error_norm.py",
                             ["density","velocity_x","velocity_y","velocity_z",
                              "total_energy","bfield_x","bfield_y","bfield_z",
                              "red"])

def run_tests(executable):
    temp = ('input/vlct/passive_advect_sound_wave/'
            'method_vlct_{}_passive_sound.in')
    call_test = EnzoEWrapper(executable,temp)

    # calls some tests
    call_test('x')
    call_test('y')
    call_test('z')

def passive_sound_l1_analyze(axis, ref_l1_norm):
    t1_dir = data_dir_template.format(axis, 0.0)
    t2_dir = data_dir_template.format(axis, 1.0)
    norm = calc_l1_norm(t1_dir,t2_dir)

    if not isclose(norm, ref_l1_norm, abs_tol = True):
        message_temp = ("L1 error for {:s}-axis aligned sound wave with "
                        "passive scalar is wrong\n{:s} {:s}")
        print(message_temp.format(axis, repr(norm), repr(ref_l1_norm)))
        return False
    return True


def analyze_tests():

    # the values are taken from runs on habanero
    # first let's run the l1-norm of each value for 2 sizes to make sure we get
    # convergence

    # originally we were looking at the slope of convergence but that allows
    # the errors to be larger

    # These errors were all obtained from parallel runs with the domain split
    # between processes
    ref_l1_error = 6.918011602177191e-08
    r = []
    r.append(passive_sound_l1_analyze('x', ref_l1_error))
    r.append(passive_sound_l1_analyze('y', ref_l1_error))
    r.append(passive_sound_l1_analyze('z', ref_l1_error))

    n_passed = np.sum(r)
    n_tests = len(r)
    print("{:d} Tests passed out of {:d} Tests.".format(n_passed,n_tests))

    return n_passed == n_tests

def cleanup():

    dir_names = ["method_vlct-1-x_passive_soundN16-0.0",
                 "method_vlct-1-x_passive_soundN16-1.0",
                 "method_vlct-1-y_passive_soundN16-0.0",
                 "method_vlct-1-y_passive_soundN16-1.0",
                 "method_vlct-1-z_passive_soundN16-0.0",
                 "method_vlct-1-z_passive_soundN16-1.0"]

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
