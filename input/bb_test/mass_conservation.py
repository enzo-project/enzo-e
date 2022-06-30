#!/bin/python

# This file defines the function test_mc which has the following behaviour:
# - takes as arguments the prefix for the names of the directories generated
#   when Enzo-E runs the bb_test problem, which contain snapshots of
#   the field data particle data at regular time intervals, a name for the
#   file to which the plot will be saved, and a value for the tolerance limit.
# - uses yt to read the data and calculates the total mass in each snapshot. It then
#   calculates the mass conservation error, and tests whether
#   the absolute value is always less than the tolerance limit.
#   If the absolute value of the error is larger than the tolerance for any
#   snapshot, then the test fails.
# - generates a plot showing the total gas mass, total sink mass, and total mass against time

# test_mc can be imported by another script (as is done in
# run_bb_test.py), or it can be executed by running this file as a script,
# with command line arguments being passed to test_mmc.
# For more information, run "python mass_conservation.py -h".

import yt
import matplotlib.pyplot as plt
import numpy as np
import io
import argparse as ap
import sys

############################################################################

def test_error(quantity_list,tolerance):

    try:
        return_val = True
        error_list = []
        initial_value = quantity_list[0]
        if np.abs(initial_value) < tolerance:
            for q in quantity_list:
                error = q - initial_value
                error_list.append(error)
                if np.abs(error) > tolerance:
                    return_val = False
        else:
            for q in quantity_list:
                error = q / initial_value - 1.0
                error_list.append(error)
                if np.abs(error) > tolerance:
                    return_val = False

        return (np.array(error_list),return_val)
    except:
        print("Encounted error while trying to calculate error in quantity")
        sys.exit(1)

def test_mc(tolerance,input_prefix,output):

    ds_pattern = f"{input_prefix}_????/{input_prefix}_????.block_list"

    time = []
    gas_mass = []
    sink_mass = []

    # Can't test conservation of quantities if there are less than 2 snapshots
    # If this has happened, something has probably gone wrong in any case, so
    # return False
    if (len(yt.DatasetSeries(ds_pattern)) < 2):
        return False

    for ds in yt.DatasetSeries(ds_pattern):
        box = ds.box(left_edge = -ds.domain_width/2.0,
                     right_edge = ds.domain_width/2.0)
        data = ds.all_data()
        gas_mass.append(data["gas","mass"].sum())
        sink_mass.append(box["sink","mass"].sum())
        time.append(ds.current_time)

    total_mass = [gm + sm for (gm,sm) in zip(gas_mass,sink_mass)]

    mass_error, mass_error_pass = test_error(total_mass,tolerance)

    print(f"Maximum mass error = {np.max(np.abs(mass_error))}")

    # replace zeros with a small number. Helpful for plotting purposes
    mass_error[mass_error == 0] = 1.0e-100

    ### make the figure

    # set matplotlib params
    params = {'axes.labelsize': 16,
              'axes.titlesize': 16,
              'font.size': 16,
              'legend.fontsize': 16,
              'xtick.labelsize': 16,
              'ytick.labelsize': 16,
              'figure.figsize' : (8,6),
              'figure.subplot.left'    : 0.15,
              'figure.subplot.right'   : 0.95  ,
              'figure.subplot.bottom'  : 0.25  ,
              'figure.subplot.top'     : 0.95  ,
              'figure.subplot.wspace'  : 0.10  ,
              'figure.subplot.hspace'  : 0.05  ,
              'lines.markersize' : 3.0,
              'lines.linewidth' : 2.0
    }

    plt.rcParams.update(params)

    fig,ax = plt.subplots()
    ax.plot(time,gas_mass/total_mass[0],label = "Gas mass")
    ax.plot(time,sink_mass/total_mass[0],label = "Sink mass")
    ax.plot(time,total_mass/total_mass[0],label = "Total mass")
    grid_line_positions = np.linspace(0.0,1.0,11,endpoint = True)
    for y in grid_line_positions:
        ax.axhline(y = y,ls = "--",color = "k",linewidth = 0.5,alpha = 0.7)

        ax.set_ylim(0.0,1.05)
    ax.set_ylabel("Mass fraction")
    ax.set_xlabel("Time (years)")
    ax.legend(loc = 0)

    fig.savefig(output)

    print(f"Saved plot to {output}")

    return mass_error_pass


if __name__ == "__main__":

    parser = ap.ArgumentParser(
    description="""
    Reads a data directory containing a series of snapshots, and
    tests whether mass is conserved.
    """
    )

    parser.add_argument(
        "-i",
        "--input_prefix",
        help="""
        The prefix for the directories containing snapshots.
        Default = \"Dir\".
        """,
        required=False,
        default="Dir",
        type=str,
    )

    parser.add_argument(
        "-o",
        "--output",
        help="Output filename. Default: ./mc.png",
        required=False,
        default="mc.png",
        type=str,
    )

    parser.add_argument(
        "-t",
        "--tolerance",
        help="""
        The error tolerance. If initial mass is greater than this value,
        this script tests if the absolute difference between the
        total mass in each snapshot and the initial mass,
        divided by the initial mass, is less than the tolerance.
        If initial mass is smaller, then just checks the absolute
        difference.
        """,
        required=True,
        type=float,
    )

    args = parser.parse_args()
    test_passed = False
    try:
        if test_mc(args.tolerance,
                   args.input_prefix,
                   args.output):
            test_passed = True

    except Exception as e:
        print("Encountered error when trying to test mass conservation:")
        print(e)
        sys.exit(2)

    if test_passed:
        print("Test passed!")
        sys.exit(0)
    else:
        print("Test failed, mass not conserved")
        sys.exit(3)
