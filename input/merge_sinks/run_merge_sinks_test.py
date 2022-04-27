#!/bin/python

# Running run_merge_sinks_test.py does the following:

# - Generates the initial conditions file called particles.dat
# - Runs Enzo-E taking merge_sinks_test.in as a parameter file, and using
#   particles.dat to make the initial conditions. This produces a set of directories
#   which contain snapshots of the particle data at regular time intervals.
# - Produces  a set of png images, which show scatter plots of the x and y coordinates
#   of particles read from the snapshot directories.
# - Analyzes snapshot data to test for mass and momentum conservation.
#   This generates a figure (mmc.png), which shows the mass conservation error,
#   the momentum conservation error, and the number of particles plotted against
#   the number of cycles (i.e., timesteps). The test is considered to have passed
#   if all the errors are less than some tolerance limit, which is 1.0e-6 in the
#   case of double precision, and 1.0e-4 for single precision.
# - Deletes particles.dat and the snapshot directories.

# run_merge_sinks_test.py takes the following arguments:

# - "--launch_cmd" which is the command used to run Enzo-E.

#   To run Enzo-E as a serial program, set --launch_cmd to /path/to/bin/enzo-e
#   To run Enzo-E as a parallel program, set --launch_cmd to
#   "/path/to/bin/charmrun +p 4 ++local /path/to/bin/enzo-e"

# - "--prec" which should be set to  "single" or "double" depending
#   on whether Enzo-E was compiled with single- or double- precision. This affects
#   the error tolerance used when testing for mass and momentum conservation, with
#   the tolerance being higher in the single-precision case.

# - "--ics_type" which can be either "stationary" or "drift". The "stationary" option
#   sets up a sphere of particles in the middle of the domain, with the particles
#   having an infall velocity but no additional drift velocity. The "drift" option
#   sets up the sphere of particles in the middle of a block, with the particles
#   having an infall velocity in addition to a uniform drift velocity.
#   "

import argparse
import os
import sys
import shutil
import subprocess
import glob

import yt

from testing_utils import testing_context
from ics import make_ics
from images import make_images
from mass_momentum_conservation import test_mmc

yt.mylog.setLevel(30) # set yt log level to "WARNING"

def run_test(executable):

    param_file = 'input/merge_sinks/merge_sinks_test.in'
    command = executable + ' ' + param_file
    subprocess.call(command,shell = True)

def analyze_test(prec):

    # set tolerance level depending on whether Enzo-E was compiled with single or
    # double precision

    if prec == "double":
        tolerance = 1.0e-6
    else:
        tolerance = 1.0e-4

    try:
        return test_mmc(tolerance,
                        input_prefix = "Dir",
                        output = "mmc.png")
    except:
        print("Encountered error when trying to test mass and momentum conservation")
        return False

def cleanup():
    dir_list = glob.glob("Dir*")
    for dir_name in dir_list:
        if os.path.isdir(dir_name):
            shutil.rmtree(dir_name)

    if os.path.isfile("particles.dat"):
        os.remove("particles.dat")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--launch_cmd', required=True,type=str)
    parser.add_argument('--prec', choices=['double', 'single'], required=True, type=str)
    parser.add_argument('--ics_type', choices=['stationary', 'drift'], required=True, type=str)
    args = parser.parse_args()

    with testing_context():

        # create initial conditions in a file called particles.dat
        if (args.ics_type == "stationary"):
            make_ics(radius = 1.0)
        else:
            make_ics(centre = [1.5,1.5,1.5],
                     drift_velocity = [2.0,2.0,2.0])

            # run the test
        run_test(args.launch_cmd)

        # make images
        make_images(input_prefix = "Dir", output_prefix = "image")

        # analyze the test
        tests_passed = analyze_test(args.prec)

        # cleanup
        cleanup()

    if tests_passed:
        sys.exit(0)
    else:
        sys.exit(3)
