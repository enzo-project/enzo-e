#!/bin/python

# runs VLCT cloud tests that explicitly check whether the implementation of the
# dual energy formalism appropriately maintain asymmetry.


import os.path
import sys
import shutil
import subprocess

import numpy as np
import yt

yt.mylog.setLevel(30) # set yt log level to "WARNING"

from testing_utils import prep_cur_dir, EnzoEWrapper

def run_tests(executable):

    temp = 'input/vlct/dual_energy_cloud/{}_cloud.in'
    wrapper = EnzoEWrapper(executable,temp)

    wrapper('hlld')
    wrapper('hllc')
    wrapper('hlle')

def slice_asym(grid, slice_ax, slice_ind, flip_across):
    """
    For a given slice through a grid, this function quantifies
    the asymmetry with respect to an axis bisecting the slice.
    """
    
    assert slice_ax in 'xyz'
    assert flip_across in 'xyz'
    assert slice_ax != flip_across
    
    if slice_ax == 'x':
        slice_arr = grid[slice_ind,:,:]
        if flip_across == 'y':
            flipped = grid[slice_ind,:,::-1]
        else:
            flipped = grid[slice_ind,::-1,:]
    elif slice_ax == 'y':
        slice_arr = grid[:,slice_ind,:]
        if flip_across == 'x':
            flipped = grid[:,slice_ind,::-1]
        else:
            flipped = grid[::-1,slice_ind,:]
    else:
        slice_arr = grid[:,:,slice_ind]
        if flip_across == 'x':
            flipped = grid[:,::-1,slice_ind]
        else:
            flipped = grid[::-1,:,slice_ind]
    
    return np.sum(np.abs((slice_arr-flipped)/slice_arr))

def check_cloud_asym(fname, name, max_asym):
    ds = yt.load(fname)
    grid = ds.covering_grid(0, ds.domain_left_edge, 
                            ds.domain_dimensions)
    density = grid["density"].in_cgs().v
    

    out = []
    # note since cloud propagates along x-axis, we flip across x twice
    for slice_ax, slice_ind, flip_ax in [('z',16,'x'),('y', 16, 'x'),
                                        ('x', 16, 'z')]:
        asym = slice_asym(density, slice_ax, slice_ind, flip_ax)
        if asym > max_asym:
            msg = ("FAILED: asymmetry of {} with respect to the central "
                   "{}-axis in the slice through the {}-axis with index={} is "
                   "{}. This exceeds the max tolerated asymmetry {}")
            print(msg.format(name,flip_ax,slice_ax,slice_ind,asym,max_asym))
            out.append(False)
        else:
            out.append(True)
    return out

def analyze_tests():
    r = []
    r += check_cloud_asym('hlld_cloud_0.0625/hlld_cloud_0.0625.block_list',
                          'hlld_cloud', 0.)
    r += check_cloud_asym('hllc_cloud_0.0625/hllc_cloud_0.0625.block_list',
                          'hllc_cloud', 0.)
    r += check_cloud_asym('hlle_cloud_0.0625/hlle_cloud_0.0625.block_list',
                          'hlle_cloud', 3.e-13)
    n_passed = np.sum(r)
    n_tests = len(r)
    print("{:d} Tests passed out of {:d} Tests.".format(n_passed,n_tests))

    return n_passed == n_tests

def cleanup():
    dir_names = ['hlld_cloud_0.0625', 'hllc_cloud_0.0625', 'hlle_cloud_0.0625']
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
