#!/bin/python

# Running run_bb_test.py does the following:

# - Runs Enzo-E ising the parameter file specified by --param_file. This produces a set
#   of directories which contain snapshots of the particle data at regular time intervals.
# - Produces  a set of png images, which show scatter plots of the x and y coordinates
#   of particles read from the snapshot directories.
# - Analyzes snapshot data to test for mass and momentum conservation.
#   This generates a figure (mmc.png), which shows the mass and momentum conservation errors.
#   The test is considered to have passed if all the errors are less than some tolerance limit,
#   which is 1.0e-12 in the case of double precision, and 1.0e-4 for single precision.
# - Deletes the snapshot directories.

# Arguments:
# --launch_cmd: the command used to run Enzo-E.
#               To run Enzo-E as a serial program, set this to `/path/to/bin/enzo-e`.
#               To run Enzo-E as a parallel program, set this to (for example)
#               `"/path/to/bin/charmrun +p4 ++local /path/to/bin/enzo-e"`
#
# --prec: Can be set to `single` or `double`, and should be the same as the precision
#         with which Enzo-E was compiled. Sets the error tolerance when testing
#         conservation of mass and momentum.
#

# If running
import argparse
import os
import sys
import shutil
import subprocess
import glob

import yt

from testing_utils import testing_context
from images import make_images
from mass_conservation import test_mc

yt.mylog.setLevel(30) # set yt log level to "WARNING"

def run_test(executable):

    command = executable + ' ' + "input/bb_test/bb.in"
    subprocess.call(command,shell = True)

def analyze_test(prec):

    # set tolerance level depending on whether Enzo-E was compiled with single or
    # double precision

    if prec == "double":
        tolerance = 1.0e-12
    else:
        tolerance = 1.0e-4

    try:
        return test_mc(tolerance,
                       input_prefix = "Dir",
                       output = "mc.png")
    except:
        print("Encountered error when trying to test mass conservation")
        return False

def cleanup():
    dir_list = glob.glob("Dir*")
    for dir_name in dir_list:
        if os.path.isdir(dir_name):
            shutil.rmtree(dir_name)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--launch_cmd', required=True,type=str)
    parser.add_argument('--prec', choices=['double', 'single'], required=True, type=str)
    args = parser.parse_args()

    with testing_context():

        # run the test
        run_test(args.launch_cmd)

        make_images(input_prefix = "Dir",
                    output_prefix = "image",
                    z_slice  = 0.0,
                    min_dens = 1.0e-18,
                    max_dens = 1.0e-14)
        
        # analyze the test
        tests_passed = analyze_test(args.prec)

        # cleanup
        cleanup()

    if tests_passed:
        sys.exit(0)
    else:
        sys.exit(3)
