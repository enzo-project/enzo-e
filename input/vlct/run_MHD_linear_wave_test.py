#!/bin/python

# runs VLCT MHD Linear Wave Tests
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


# this executes things in standalone mode
l1_norm_calc_template = ("python tools/l1_error_norm.py sim {:s} {:s} -n {:d}"
                         " -f density,velocity_x,velocity_y,velocity_z,"
                         "pressure,bfield_x,bfield_y,bfield_z")
data_dir_template = "method_vlct-{:d}-{:s}N{:d}_{:.1f}"

def run_tests(executable):
    temp = 'input/vlct/MHD_linear_wave/method_vlct_{}N{:d}{}.in'
    call_test = EnzoEWrapper(executable,temp)

    call_test("fast",16,'')
    call_test("fast",32,'')
    call_test("fast",32,'_left')

    call_test("alfven",16,'')
    call_test("alfven",32,'')

    call_test("slow",16,'')
    call_test("slow",32,'')
    
    call_test("entropy",16,'')
    call_test("entropy",32,'')

def format_pair(template, nblocks, res, wave_name, t1=0., t2=0., direction=''):
    # formats a pair of directory templates
    return (template.format(nblocks=nblocks, res=res, wave_name=wave_name,
                            direction=direction, time=t1),
            template.format(nblocks=nblocks, res=res, wave_name=wave_name,
                            direction=direction, time=t2))

def analyze_linwave(ref_l1_norm, wave_name, nblocks, res, target_template,
                    name_template, l1_functor, final_times, verbose = False):
    """
    Both name_template and target_template are used to format strings based on 
    other arguments with the `format` method. They should each employ named 
    placeholders and expect `format` to be passed the following kwargs:
        `direction`, `wave_name`, `res`, and `nblocks`.
    Note that `direction` is set to ''. Additionally, target_template should 
    also expect a `time` kwarg
    """
    common = dict(nblocks=nblocks, res=res, wave_name=wave_name, direction='')
    t1_dir, t2_dir = format_pair(target_template, t1 = 0.,
                                 t2 = final_times[wave_name], **common)
    test_case_name = '{wave_name} wave N={res}'.format(**common)
    functor_kwargs = dict(dir_name2 = t2_dir, res = res)

    return standard_analyze(ref_l1_norm, l1_functor, target_path = t1_dir,
                            test_case_name = test_case_name, verbose = verbose,
                            functor_kwargs = functor_kwargs)


# part of the following should probably be folded into standard_analyze
def identical_l1_error_linwave(nblocks, wave_name, res, l1_functor, template,
                               final_time, nblocks2 = None, prec = None,
                               verbose = False):
    """
    Compares linear wave l1 norms that should be identical.

    There are 2 main cases where this is expected:
      1. Left propagation vs. right propation
      2. Changing the number of blocks over which the simulation is distibuted.

    Parameters
    ----------
    nblocks: int
        Number of blocks over which the primary simulation's grid is divided
    wave_name: string
        The wave type - "fast", "slow", "alfven", "entropy"
    res: int
        The resolution of the simulation. For linear waves, a value of N 
        implies (2N,N,N) for (x,y,z).
    l1_functor : instance of CalcSimL1Norm
        This functor wraps the call to the l1_error_norm script.
    template : str
        Formatting this template string, should provide paths to directories 
        holding the simulation data that are used to compute the norms of L1 
        error vectors. The templates should used named placeholders for 
        formatting. This function provides the following placeholders:
            `wave_name`, `res`, `nblocks`, `direction`, and `time`
        The first 3 placeholders will be substituted with the appropriate 
        values passed by this function. The `direction` placeholder will only 
        be substituted with `''` or `'-left'`. Finally, the `time` placeholder
        will be substituted by `0.0` and the value returned by 
        `final_time[wave_type]`
    final_time : float
        The time at which the simulation of the wave ends
    nblocks2: int, optional
        Number of blocks over which the secondary simulation's grid is divided.
        If this value is not specified, then the number of blocks are assumed
        to be the same. In this case, the wave in the primary (secondary) 
        simulation is assumed to propagate rightwards (leftwards)
    prec: int, optional
        Number of digits after the decimal that must match. Athena++ 
        traditionally just compares 6. If this is passed None (default), then
        all digits are compared.
    verbose : bool, optional
        If true, prints note about tests that have passed
    """

    fmt_kwargs = dict(nblocks = nblocks, res = res, wave_name = wave_name,
                      t1 = 0., t2 = final_time, direction = '')
    t0_dir,tf_dir = format_pair(template, **fmt_kwargs)
    ref_norm = l1_functor(t0_dir,tf_dir,res)

    if nblocks2 is None:
        fmt_kwargs['direction'] = '-left'
        insert = 'right & left propagating {wave} waves'
    else:
        fmt_kwargs['nblocks'] = nblocks2
        insert = "{wave} waves split over {nblocks} and {nblocks2} blocks"

    t0_dir,tf_dir = format_pair(template, **fmt_kwargs)
    comp_norm = l1_functor(t0_dir,tf_dir,res)

    msg = None
    if prec is not None:
        ref_norm  = '{:.{prec}e}'.format(ref_norm,  prec=prec)
        comp_norm = '{:.{prec}e}'.format(comp_norm, prec=prec)
    passed = (ref_norm == comp_norm)
    if not passed:
        msg = "FAILED: L1 norms of " + insert + " don't match\n{ref} {comp}"
    elif verbose:
        msg = "PASSED: L1 norms of " + insert + " are both {ref}"

    if msg is not None:
        print(msg.format(nblocks=nblocks, nblocks2=nblocks2, wave = wave_name,
                         ref=repr(ref_norm), comp=repr(comp_norm)))
    return passed

def analyze_tests():
    # define the functor for evaluating the norm of the L1 error vector
    l1_func = CalcSimL1Norm("tools/l1_error_norm.py",
                            ["density","velocity_x","velocity_y","velocity_z",
                             "pressure","bfield_x","bfield_y","bfield_z"])

    # define the template for the directory holding the simulation data
    template = "method_vlct-{nblocks}-{wave_name}N{res}{direction}_{time:.1f}"

    final_times = {"fast" : 0.5, "alfven" : 1.0, "entropy" : 1.0, "slow" : 2.0}

    err_compare = partial(analyze_linwave, target_template = template,
                          name_template = '{wave_name} wave N={res}',
                          l1_functor = l1_func, final_times = final_times,
                          verbose = False)

    r = []
    # first let's run the l1-norm of each value for 2 sizes to make sure we get
    # convergence (originally we were looking at the slope of convergence but
    # that allows the errors to be larger)

    r.append(err_compare(1.6388526155394664e-07, 'fast',    1, 16))
    r.append(err_compare(3.302538226654406e-08, 'fast',    1, 32))

    r.append(err_compare(1.927245356389947e-07, 'alfven',  1, 16))
    r.append(err_compare(3.005870212811956e-08,  'alfven',  1, 32))

    r.append(err_compare(2.2373810027584788e-07, 'slow',    1, 16))
    r.append(err_compare(4.437025228314115e-08, 'slow',    1, 32))

    r.append(err_compare(1.0021263485338544e-07, 'entropy', 1, 16))
    r.append(err_compare(2.9194839706868883e-08,  'entropy', 1, 32))

    # Check error between left and right propagating waves
    print("The following test has never passed in the history of this "
          "implementation")
    r.append(identical_l1_error_linwave(1, 'fast', 32, l1_func, template, 0.5,
                                        verbose = True))

    n_passed = np.sum(r)
    n_tests = len(r)
    success = (n_passed == n_tests)
    print("{:d} Tests passed out of {:d} Tests.".format(n_passed,n_tests))
    if (np.sum(r[1:]) == (n_tests-1)):
        print("All tests that the this VL+CT implementation has ever passed, "
              "have suceeded")
        success = True
    return success

def cleanup():

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
