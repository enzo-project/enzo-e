#!/bin/python

# runs VLCT sod Shock Tube Tests that explicitly try to check whether the
# dual energy formalsim is operating correctly. A Galilean transformation has
# been performed on the problem so that in the frame of the mesh, the shock
# tube setup is moving at a supersonic speed (Mach Number ~ 10)
# - This primarily tests that implementation of the dual energy formalism is
#   correct. These tests caught an indexing error
#
# Notes:
# - This test doesn't rigorously test the use of the dual energy formalism in
#   a case where it is essential. This would require a background velocity
#   with a larger mach number larger Mach number.
# - This test is somewhat slow. Consequentially, I decreased the resolution
#   from cell-widths of 1/256 (used in the initial conception of the test) to
#   1/128. After this change, 4 cores on my desktop to can complete the test in
#   ~1.5 min.
#
# In the future:
# - This should probably be transitioned to testing framework once it's created
# - This script expects to be called from the root level of the repository
#   OR at the same level where its defined


import os
import os.path
import sys
import numpy as np
import shutil
import subprocess
import math

from testing_utils import prep_cur_dir, EnzoEWrapper, CalcTableL1Norm
from run_shock_tube_test import analyze_shock

def run_tests(executable):

    temp = 'input/vlct/dual_energy_shock_tube/method_vlct_sod_{}_de_M10{}.in'
    wrapper = EnzoEWrapper(executable,temp)

    wrapper('x','')
    wrapper('x','_reverse')
    wrapper('y','')
    wrapper('z','')

def analyze_tests():

    # set up some basic stuff before analyzing results:
    ref_table_path = ('input/vlct/dual_energy_shock_tube/'
                      'sod_shock_tube_t0.25_res128.csv')
    verbose_l1_calc = False # Print out command line call
    verbose_analyze = False # Print out summaries of passed tests
    l1_functor = CalcTableL1Norm("tools/l1_error_norm.py",
                                 default_ref_table = ref_table_path,
                                 verbose = verbose_l1_calc)
    kwargs = dict(
        target_template = 'method_vlct-1-sod_{}_de_M10_0.25/',
        name_template = 'dual_energy_sod_{}', l1_functor = l1_functor,
        bkg_vel = [11.875,0,0], offset = 2.96875, verbose = verbose_analyze
        )

    r = []
    # check the average L1-Norm along the active axis (averaged over multiple
    # slices along the active axis). Also check that standard deviations
    # directions perpendicular to the axis of evolution are zero

    expected_err = 0.02605861216339738

    r.append(analyze_shock(expected_err, 'z', **kwargs))
    r.append(analyze_shock(0.0, 'z', std_dev = True, **kwargs))
    r.append(analyze_shock(expected_err, 'y', **kwargs))
    r.append(analyze_shock(0.0, 'y', std_dev = True, **kwargs))
    r.append(analyze_shock(expected_err, 'x', **kwargs))
    r.append(analyze_shock(0.0, 'x', std_dev = True, **kwargs))

    r.append(analyze_shock(
        expected_err, 'x', name_template = 'dual_energy_sod_{}_left',
        target_template = 'method_vlct-1-sod_x_de_M10_reverse_0.25/',
        reverse = True, l1_functor = l1_functor, bkg_vel = (-11.875,0,0),
        offset = -2.96875, verbose = verbose_analyze))

    n_passed = np.sum(r)
    n_tests = len(r)
    print("{:d} Tests passed out of {:d} Tests.".format(n_passed,n_tests))

    return n_passed == n_tests

def cleanup():

    dir_names = ["method_vlct-1-sod_z_de_M10_0.25",
                 "method_vlct-1-sod_y_de_M10_0.25",
                 "method_vlct-1-sod_x_de_M10_0.25",
                 "method_vlct-1-sod_x_de_M10_reverse_0.25"]
    for dir_name in dir_names:
        if os.path.isdir(dir_name):
            shutil.rmtree(dir_name)

if __name__ == '__main__':

    if len(sys.argv) == 1:
        nproc = 1
    elif len(sys.argv) == 2:
        try:
            nproc = int(sys.argv[1])
        except:
            raise ValueError(
                "{} can take an optional positive integer ".format(sys.argv[0])
                + "argument for specifying the use of multiple processes")
    else:
        raise ValueError(
            "{} can take up to 1 argument to specify the ".format(sys.argv[0])
            + "number of processes. {:d} arguments ".format(len(sys.argv)-1)
            + "were provided")

    if nproc <= 0:
        raise ValueError("number of processes must be a positive integer")
    elif nproc == 1:
        executable = 'bin/enzo-p'
    else:
        charm_args = os.environ.get('CHARM_ARGS')
        if charm_args is None:
            executable_template = 'charmrun +p{:d} bin/enzo-p'
        else:
            executable_template \
                = ' '.join(['charmrun', charm_args, '+p{:d}', 'bin/enzo-p'])
        executable = executable_template.format(nproc)

    # this script can either be called from the base repository or from
    # the subdirectory: input/vlct
    prep_cur_dir('bin/enzo-p')

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
