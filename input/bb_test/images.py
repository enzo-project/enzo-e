#!/bin/python

# This file defines a function called make_images, which
# takes as input the directories generated when Enzo-E runs the test
# problem, which contain snapshots of the particle data and field data
# at regular time intervals.
# It then uses yt make slice images of the density field, also showing the
# position of the sink particle. A series of image files are created, with
# each image corresponding to a particular snapshot. These images are useful
# to check whether the
# initial conditions were set up correctly and if Enzo-E ran as expected.

# make_images can be imported by another script (as is done in
# run_bb_test.py), or it can be executed by running this file
# as a script, with command line arguments being passed to make_images.
# For more information, run "python images.py -h".

import yt
import argparse as ap

def make_images(input_prefix,output_prefix,z_slice,min_dens,max_dens):

    ds_pattern = f"{input_prefix}_????/{input_prefix}_????.block_list"
    center = [0.0,0.0,z_slice]
    for ds in yt.DatasetSeries(ds_pattern):

        slc = yt.SlicePlot(ds,"z",("gas","density"),
                           center = center,
                           width = ds.domain_width[0])

        box = ds.box(left_edge = -ds.domain_width/2.0,
                     right_edge = ds.domain_width/2.0)
        sink_positions = [s for s in zip(box["sink","x"],
                                         box["sink","y"],
                                         box["sink","z"])]
        for s in sink_positions:
            slc.annotate_marker(s,
                                coord_system = "data",
                                plot_args = {'color':"k"})

        slc.annotate_grids(linewidth = 0.5,alpha = 1.0)
        slc.set_zlim(("gas","density"),min_dens,max_dens)

        current_cycle = ds.parameters["current_cycle"]
        filename = f"{output_prefix}_{current_cycle:03}.png"
        slc.save(filename)

if __name__ == "__main__":

    parser = ap.ArgumentParser(
        description="""
        Reads a data directory containing a series of snapshots of the
        Shu Collapse problem and makes an image showing a slice of the
        density field and the x and y position of the sink particle,
        making a separate image for each snapshot.
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

    parser.add_argument(
        "-z",
        "--z_slice",
        help="""
        The z-coordinate of the density slice.
        """,
        required=True,
        type=float,
    )

    parser.add_argument(
        "--min_dens",
        help="""
        The minimum density shown in the images.
        """,
        required=True,
        type=float,
    )
    parser.add_argument(
        "--max_dens",
        help="""
        The maximum density shown in the images.
        """,
        required=True,
        type=float,
    )

    args = parser.parse_args()
    make_images(args.input_prefix,args.output_prefix,args.z_slice)
