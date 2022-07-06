#!/bin/python

# This file defines a function called make_images, which
# takes as input the directories generated when Enzo-E runs the test
# problem, which contain snapshots of the particle data at regular time intervals.
# It then uses yt to read the data and then uses matplotlib to make scatter
# plots of the x and y coordinates of the particles, outputting a series of
# image files, with each image corresponding to a particular snapshot. These
# images are useful to check whether the initial conditions were set up
# correctly and if Enzo-E ran as expected.

# make_images can be imported by another script (as is done in
# run_merge_sinks_test.py), or it can be executed by running this file as a script,
# with command line arguments being passed to make_images. 
# For more information, run "python images.py -h".

import yt
import matplotlib.pyplot as plt
import argparse as ap
import matplotlib
matplotlib.use("Agg")

def make_images(input_prefix,output_prefix):
    
    ds_pattern = f"{input_prefix}_????/{input_prefix}_????.block_list"
    for ds in yt.DatasetSeries(ds_pattern):
        
        box = ds.box(left_edge = -ds.domain_width/2.0,
                     right_edge = ds.domain_width/2.0)
        x = box["sink","x"]
        y = box["sink","y"]
        z = box["sink","z"]
        fig,ax = plt.subplots()
        ax.plot(x,y,marker = "x",linestyle = "None")
        ax.set_xlim(-ds.domain_width[0]/2.0,ds.domain_width[0]/2.0)
        ax.set_ylim(-ds.domain_width[0]/2.0,ds.domain_width[0]/2.0)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        current_cycle = ds.parameters["current_cycle"]
        filename = f"{output_prefix}_{current_cycle:03}.png"
        fig.savefig(filename)
        plt.close(fig)

if __name__ == "__main__":
    
    parser = ap.ArgumentParser(
        description="""
        Reads a data directory containing a series of snapshots of the MergeSinksTest
        problem and plots the x and y positions of the sink particles,
        making a separate image for each snapshot
        """
    )

    parser.add_argument(
        "-i",
        "--input_prefix",
        help="The prefix for the directories containing snapshots",
        required=True,
        type=str,
    )

    parser.add_argument(
        "-o",
        "--output_prefix",
        help="""
        Prefix for image filenames.
        """,
        required=True,
        type=str,
    )

    args = parser.parse_args()
    make_images(args.input_prefix,args.output_prefix)
