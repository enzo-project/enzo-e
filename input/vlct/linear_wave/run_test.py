#!/bin/python

# runs VLCT Linear Wave Tests
# - this should probably be transitioned to the SConscript file in enzo-e/test.
#   However, I'm unclear on how to do that
# - This script expects to be called from the root directory of the repository
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
#
# Optional Test:
#   At one point I was slightly concerned about differences between waves
#   evolved serially and in parallel (However after testing it - this doesn't
#   seem to be a problem). Because it's difficult to launch a parallel
#   simulation from this script on a personal machine - this test optional
#
#   To run this test: prior to calling this script, run the parallel simulation
#   by calling the following command from the root directory of the repository:
#    $ charmrun ++local bin/enzo-p input/vlct/linear_wave/method_vlct_fastN32_parallel.in


import os
import os.path
import numpy as np
import shutil
import subprocess
import math


# this executes things in standalone mode
_executable = 'bin/enzo-p'
l1_norm_calc_template = "python tools/l1_error_norm.py {:s} {:s} -n {:d}"
data_dir_template = "method_vlct-{:d}-{:s}N{:d}_{:.1f}"
final_times = {"fast" : 0.5, "alfven" : 1.0, "entropy" : 1.0, "slow" : 2.0}

def call_test(wave,res,left=False):

    input_file = 'input/vlct/linear_wave/method_vlct_{}N{:d}.in'.format(wave,
                                                                        res)
    if left:
        input_file = ('input/vlct/linear_wave/' +
                      'method_vlct_{}N{:d}_left.in'.format(wave,res))
    command = _executable + ' ' + input_file
    subprocess.call(command,shell=True)

def run_tests():

    # can't actually run the serial vs parallel test from this script (not sure
    # how to get around request for password)
    call_test("fast",16)
    call_test("fast",32)
    call_test("fast",32,left=True)

    call_test("alfven",16)
    call_test("alfven",32)

    call_test("slow",16)
    call_test("slow",32)
    
    call_test("entropy",16)
    call_test("entropy",32)

def calc_l1_norm(dir1,dir2,res):
    command = l1_norm_calc_template.format(dir1, dir2, res)
    return float(subprocess.check_output(command,shell=True))

def typical_l1_norm(num_blocks, wave_name, res, time_1, time_2):
    t1_dir = data_dir_template.format(num_blocks, wave_name, res, time_1)
    t2_dir = data_dir_template.format(num_blocks, wave_name, res, time_2)
    return calc_l1_norm(t1_dir,t2_dir,res)

def isclose(a,b,abs_tol = False):
    # just wraps numpy's isclose function
    # this is how much we allow things to be wrong (a little arbitrary)
    err = 1.e-13
    # it makes sense to allow for absolute errors in l1 norms since the
    # floating point differences will enter into the values that are being
    # subtracted
    if abs_tol:
        # slow wave requires slightly bigger tolerance (since it includes more
        # timesteps which allows the floating point errors to grow more)
        atol = 1.e-14
    else:
        atol = 0
    return np.isclose(a,b,rtol=err,atol=atol)

def standard_l1_analyze(num_blocks,wave_name, res, ref_l1_norm,):
    # just going to factor out some code for some standard l1 Norm tests

    norm = typical_l1_norm(num_blocks, wave_name, res, 0.0,
                           final_times[wave_name])

    if not isclose(norm, ref_l1_norm, abs_tol = True):
        message_temp = "L1 error of {:s} wave N={:d} isn't correct\n{:s} {:s}"
        print(message_temp.format(wave_name, res, repr(norm),
                                  repr(ref_l1_norm)))
        return False
    return True

def identical_l1_error_comp(num_blocks,wave_name,res, num_blocks2 = None):
    """
    Compares l1 norms that should be identical.

    There are 2 main cases where this is expected:
      1. Left propagation vs. right propation
      2. Changing the number of blocks over which the simulation is distibuted.

    Parameters
    ----------
    num_blocks: int
        Number of blocks over which the domain of the primary simulation is 
        divided
    wave_name: string
        The wave type - "fast", "slow", "alfven", "entropy"
    res: int
        The resolution of the simulation. For linear waves, a value of N 
        implies (2N,N,N) for (x,y,z).
    num_blocks2: int, optional
        The number of blocks over which the domain of the secondary simulation 
        is divided. If this value is not specified, then the number of blocks 
        are assumed to be the same. In this case, the wave in the primary 
        (secondary) simulation is assumed to propagate rightwards (leftwards)
    """

    tf = final_times[wave_name]
    ref_norm = typical_l1_norm(num_blocks, wave_name, res, 0.0, tf)

    dir_template = "method_vlct-{:d}-fastN{:d}{:s}_{:.1f}"
    if num_blocks2 is None:
        t0_dir = dir_template.format(num_blocks,res,'-left',0.0)
        tf_dir = dir_template.format(num_blocks,res,'-left',tf)
    else:
        t0_dir = dir_template.format(num_blocks2,res,'',0.0)
        tf_dir = dir_template.format(num_blocks2,res,'',tf)
    comp_norm = calc_l1_norm(t0_dir,tf_dir,res)

    if ref_norm == comp_norm:
        return True

    if num_blocks2 is None:
        message_temp = ("L1 norms of right & left propagating {:s} waves don't "
                        "match\n{:s}, {:s}")
        message = message_temp.format(wave_name, repr(ref_norm),
                                      repr(comp_norm))
    else:
        message_temp = ("L1 norms of {:s} wave split over {:d} and {:d} blocks "
                        "don't match\n{:s}, {:s}")
        message = message_temp.format(wave_name, num_blocks, num_blocks2,
                                      repr(ref_norm), repr(comp_norm))
    print(message)
    return False

def analyze_tests():

    # the values are taken from runs on habanero
    # first let's run the l1-norm of each value for 2 sizes to make sure we get
    # convergence

    # originally we were looking at the slope of convergence but that allows
    # the errors to be larger

    # These errors were all obtained from parallel runs with the domain split
    # between processes
    r = []
    r.append(standard_l1_analyze(1,"fast", 16, 1.6388526155538374e-07))
    r.append(standard_l1_analyze(1,"fast", 32, 3.3025382238175106e-08))

    r.append(standard_l1_analyze(1,"alfven", 16, 1.9272453542097202e-07))
    r.append(standard_l1_analyze(1,"alfven", 32, 3.005870203721324e-08))
    
    r.append(standard_l1_analyze(1,"slow", 16, 2.2373810031907528e-07))
    r.append(standard_l1_analyze(1,"slow", 32, 4.437024386763274e-08))

    r.append(standard_l1_analyze(1,"entropy", 16, 1.0021263478396943e-07))
    r.append(standard_l1_analyze(1,"entropy", 32, 2.9194839306558322e-08))

    # Check error between left and right propagating waves
    r.append(identical_l1_error_comp(1, "fast", 32))

    n_passed = np.sum(r)
    n_tests = len(r)
    print("{:d} Tests passed out of {:d} Tests.".format(n_passed,n_tests))

    # Optionally, check error between serial and linear propagation
    # I can't quite figure out how to start the parallel test from this script
    # for an SMP installation without requesting a password

    if (os.path.isdir("method_vlct-8-fastN32_0.0") and
        os.path.isdir("method_vlct-8-fastN32_0.5")):
        print("Optional serial vs. parallel test.")
        if identical_l1_error_comp(1, "fast", 32, 8):
            print("Passed")

def cleanup():
    #"method_vlct-1-fastN32_cycle55", #serial vs. parallel
    #"method_vlct-8-fastN32_cycle55", #serial vs. parallel

    dir_names = ["method_vlct-1-fastN16_0.0",
                 "method_vlct-1-fastN16_0.5",
                 "method_vlct-1-fastN32_0.0",
                 "method_vlct-1-fastN32_0.5",
                 "method_vlct-1-fastN32-left_0.0",
                 "method_vlct-1-fastN32-left_0.5",
                 "method_vlct-1-alfvenN16_0.0",
                 "method_vlct-1-alfvenN16_1.0",
                 "method_vlct-1-alfvenN32_0.0",
                 "method_vlct-1-alfvenN32_1.0",
                 "method_vlct-1-slowN16_0.0",
                 "method_vlct-1-slowN16_2.0",
                 "method_vlct-1-slowN32_0.0",
                 "method_vlct-1-slowN32_2.0",
                 "method_vlct-1-entropyN16_0.0",
                 "method_vlct-1-entropyN16_1.0",
                 "method_vlct-1-entropyN32_0.0",
                 "method_vlct-1-entropyN32_1.0"]

    for dir_name in dir_names:
        if os.path.isdir(dir_name):
            shutil.rmtree(dir_name)

if __name__ == '__main__':

    # this script can either be called from the base repository or from
    # the subdirectory: input/vlct/linear_wave
    cwd = os.getcwd()
    if cwd[-22:] == "input/vlct/linear_wave":
        os.chdir("../../../")
        # the following is just for a possible error message
        orig_exec_path = os.path.join(cwd,"../../../",_executable)
    else:
        # the following is just for a possible error message
        orig_exec_path = os.path.join(cwd,_executable)

    if not os.path.isfile(_executable):
        raise RuntimeError("Can't locate the executable: " + orig_exec_path)

    # run the tests
    run_tests()

    # analyze the tests
    analyze_tests()

    # cleanup the tests
    cleanup()
