#!/bin/python

# This file defines the function test_mmc which has the following behaviour:
# - takes as the prefix for the names of the directories generated when Enzo-E runs
#   the accretion_test problem, which contain snapshots of the
#   particle data at regular time intervals, a name for the image file which will be
#   generated, and a value for the tolerance limit
# - uses yt to read the data and calculates the total mass, total x/y/z-momentum
#   and total number of particles in each snapshot. It then calculates the
#   conservation error for these quantities, and tests whether the absolute value
#   is less than the tolerance limit. If the absolute value of the error is larger
#   than the tolerance for any snapshot, then the test fails.
# - generates a figure with three panels, showing the mass conservation error,
#   momentum conservation error, and number of particles plotted against the number
#   of cycles (i.e. timesteps).  This figure is useful for checking if particles are
#   indeed merging (shown by decrease in particle number with time) and whether
#   mass and momentum are being properly conserved.

# test_mmc can be imported by another script (as is done in run_accretion_test.py), 
# or it can be executed by running this file as a script, with command line arguments 
# being passed to test_mmc.
# For more information, run "python mass_momentum_conservation.py -h".

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

def test_mmc(tolerance,input_prefix,output):
    
    ds_pattern = f"{input_prefix}_????/{input_prefix}_????.block_list"

    cycle = []
    gas_mass = []
    gas_px   = []
    gas_py   = []
    gas_pz   = []
    star_pz   = []
    star_mass = []
    star_px   = []
    star_py   = []
    star_pz   = []
    
    for ds in yt.DatasetSeries(ds_pattern):
        box = ds.box(left_edge = -ds.domain_width/2.0,
                     right_edge = ds.domain_width/2.0)
        data = ds.all_data()
        gas_mass.append(data["gas","mass"].sum())
        gas_px.append((data["gas","mass"] * data["gas","velocity_x"]).sum())
        gas_py.append((data["gas","mass"] * data["gas","velocity_y"]).sum())
        gas_pz.append((data["gas","mass"] * data["gas","velocity_z"]).sum())
        star_mass.append(box["star","mass"].sum())
        star_px.append((box["star","mass"] * box["star","vx"]).sum())
        star_py.append((box["star","mass"] * box["star","vy"]).sum())
        star_pz.append((box["star","mass"] * box["star","vz"]).sum())
        cycle.append(ds["current_cycle"])

    total_mass = [gm + sm for (gm,sm) in zip(gas_mass,star_mass)]
    total_px = [gpx + spx for (gpx,spx) in zip(gas_px,star_px)]
    total_py = [gpy + spy for (gpy,spy) in zip(gas_py,star_py)]
    total_pz = [gpz + spz for (gpz,spz) in zip(gas_pz,star_pz)]

    mass_error, mass_error_pass = test_error(total_mass,tolerance)
    px_error, px_error_pass = test_error(total_px,tolerance)
    py_error, py_error_pass = test_error(total_py,tolerance)
    pz_error, pz_error_pass = test_error(total_pz,tolerance)

    print(f"Maximum mass error = {np.max(np.abs(mass_error))}")
    print(f"Maximum x momentum error = {np.max(np.abs(px_error))}")
    print(f"Maximum y momentum error = {np.max(np.abs(py_error))}")
    print(f"Maximum z momentum error = {np.max(np.abs(pz_error))}")

    ### make the figure
    
    # set matplotlib params
    params = {'axes.labelsize': 16,
              'axes.titlesize': 16,
              'font.size': 16,
              'legend.fontsize': 16,
              'xtick.labelsize': 16,
              'ytick.labelsize': 16,
              'figure.figsize' : (7,9),
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

    err_lim_plot = 10.0 * tolerance
    mass_label = f"Total mass error / {err_lim_plot:.2g}"
    momentum_label = f"Total momentum error / {err_lim_plot:.2g}"
    fig,ax = plt.subplots(nrows = 2,sharex = True)
    grid_line_positions = np.linspace(-1.0,1.0,11,endpoint = True)
    
    for y in grid_line_positions:
        ax[0].axhline(y = y,ls = "--",color = "k",linewidth = 0.5,alpha = 0.7)
        ax[1].axhline(y = y,ls = "--",color = "k",linewidth = 0.5,alpha = 0.7)

    
    ax[0].plot(cycle, mass_error / err_lim_plot)
    ax[0].axhline(y = 0.1, ls = "--", color = "r", alpha = 0.7)
    ax[0].axhline(y = -0.1, ls = "--", color = "r", alpha = 0.7)
    ax[0].set_ylim(-1.0,1.0)
    ax[0].set_ylabel(mass_label)
    ax[1].plot(cycle,px_error / err_lim_plot,label = "x-momentum")
    ax[1].plot(cycle,py_error / err_lim_plot,label = "y-momentum",ls = "--")
    ax[1].plot(cycle,pz_error / err_lim_plot,label = "z-momentum",ls = ":")
    ax[1].axhline(y = 0.1, ls = "--", color = "r", alpha = 0.7)
    ax[1].axhline(y = -0.1, ls = "--", color = "r", alpha = 0.7,label = "Tolerance")
    ax[1].set_ylim(-1.0,1.0)
    ax[1].set_ylabel(momentum_label)
    ax[1].legend(loc = "upper left",bbox_to_anchor = (0.0,-0.2),ncol = 2)
    ax[1].set_xlabel("Cycle")
    
    fig.savefig(output)
    plt.close(fig)

    print(f"Saved plot to {output}")
    
    return mass_error_pass and px_error_pass and py_error_pass and pz_error_pass

if __name__ == "__main__":

    parser = ap.ArgumentParser(
    description="""
    Reads a data directory containing a series of snapshots, and
    makes plots showing mass and momentum conservation
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
        "--output",
        help="Output filename. Default: ./mmc.png",
        required=False,
        default="mmc.png",
        type=str,
    )

    parser.add_argument(
        "-t",
        "--tolerance",
        help="""
        The error tolerance. If initial mass / momentum is greater than this value,
        this script tests if the absolute difference between the 
        total mass / momentum in each snapshot and the initial mass / momentum, 
        divided by the initial mass / momentum, is less than the tolerance.
        If initial mass / momentum is smaller, then just check the absolute 
        difference.
        """,
        required=True,
        type=float,
    )

    args = parser.parse_args()
    test_passed = False
    try:
        if test_mmc(args.tolerance,
                    args.input_prefix,
                    args.output):
            test_passed = True

    except:
        print("Encountered error when trying to test mass and momentum conservation")
        sys.exit(2)

    if test_passed:
        print("Test passed!")
        sys.exit(0)
    else:
        print("Test failed, mass and momentum not conserved")
        sys.exit(3)
        

    
