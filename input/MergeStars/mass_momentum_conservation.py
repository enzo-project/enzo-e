# mass_momentum_conservation.py takes as input the directories
# generated when Enzo-E runs the test problem, which contain snapshots of the
# particle data at regular time intervals. The script uses yt to read the
# data and calculates the total mass, total x/y/z-momentum and total number
# of particles in each snapshot. It then calculates the conservation error for these
# quantities, and tests whether the absolute value is less than some tolerance,
# which is set to 1e-6 by default. If the absolute value of the error is larger
# than the tolerance for any snapshot, then the script returns an error code.

# Finally, the script generates a figure
# with three panels, showing the mass conservation error, momentum conservation
# error, and number of
# particles plotted against the number of cycles (i.e. timesteps). The mass
# / momentum error at a given cycle is defined as the difference between the total
# mass / momentum at a given cycle and the initial total mass / momentum,
# divided by the initial total mass / momentum. This figure is useful for
# checking if particles are indeed merging (shown by decrease in particle number
# with time) and whether mass and momentum are being properly conserved.
# This script is designed to be run
# with MPI as "mpirun -np X python mass_momentum_conservation.py", with X
# the number of processes.
# For further details run "python mass_momentum_conservation.py -h".

import yt
import matplotlib.pyplot as plt
import numpy as np
import io
import argparse as ap
import sys

# matplotlib params
params = {'axes.labelsize': 16,
          'axes.titlesize': 16,
          'font.size': 16,
          'legend.fontsize': 16,
          'xtick.labelsize': 16,
          'ytick.labelsize': 16,
          'text.usetex': True,
          'figure.figsize' : (7,12),
          'figure.subplot.left'    : 0.15,
          'figure.subplot.right'   : 0.95  ,
          'figure.subplot.bottom'  : 0.25  ,
          'figure.subplot.top'     : 0.95  ,
          'figure.subplot.wspace'  : 0.10  ,
          'figure.subplot.hspace'  : 0.05  ,
          'lines.markersize' : 3.0,
          'lines.linewidth' : 2.0
          #'text.latex.unicode': True
      }

plt.rcParams.update(params)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Times']})

############################################################################

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
    this script tests if the absolute difference between the total mass / momentum
    in each snapshot and the initial mass / momentum, divided by the initial mass /
    momentum, is less than the tolerance. If initial mass / momentum is smaller,
    then just check the absolute difference.
    """,
    required=False,
    default=1.0e-6,
    type=float,
)

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

        return (error_list,return_val)
    except:
        print("Encounted error while trying to calculate error in quantity")
        sys.exit(2)
 
if __name__ == "__main__":
    
    args = vars(parser.parse_args())
    input_prefix = args["input_prefix"]
    yt.enable_parallelism()
    ds_pattern = f"{input_prefix}_????/{input_prefix}_????.block_list"
    ts =  yt.DatasetSeries(ds_pattern)

    storage = {}
    for sto,ds in ts.piter(storage = storage):
        box = ds.box(left_edge = -ds.domain_width/2.0,right_edge = ds.domain_width/2.0)
        mass = box["star","mass"].sum()
        px = (box["star","mass"] * box["star","vx"]).sum()
        py = (box["star","mass"] * box["star","vy"]).sum()
        pz = (box["star","mass"] * box["star","vz"]).sum()
        n_p = len(box["star","mass"])
        sto.result = ds.parameters["current_cycle"],mass,px,py,pz,n_p

    time_list = []
    mass_list = []
    px_list = []
    py_list = []
    pz_list = []
    n_p_list = []
    if yt.is_root():
        for i, (t,mass,px,py,pz,n_p) in sorted(storage.items()):
            time_list.append(t)
            mass_list.append(mass)
            px_list.append(px)
            py_list.append(py)
            pz_list.append(pz)
            n_p_list.append(n_p)

        mass_error, mass_error_pass = test_error(mass_list,args["tolerance"])
        px_error, px_error_pass = test_error(px_list,args["tolerance"])
        py_error, py_error_pass = test_error(py_list,args["tolerance"])
        pz_error, pz_error_pass = test_error(pz_list,args["tolerance"])
        
        fig,ax = plt.subplots(nrows = 3,sharex = True)
        ax[0].plot(time_list,mass_error)
        grid_line_positions = np.linspace(-0.1,0.1,11,endpoint = True)
        for y in grid_line_positions:
            ax[0].axhline(y = y,ls = "--",color = "k",linewidth = 0.5,alpha = 0.7)
            ax[1].axhline(y = y,ls = "--",color = "k",linewidth = 0.5,alpha = 0.7)

        ax[0].set_ylim(-0.1,0.1)
        ax[0].set_ylabel("Total Mass Error")

        ax[1].plot(time_list,px_error,label = r"$p_x$")
        ax[1].plot(time_list,py_error,label = r"$p_y$",ls = "--")
        ax[1].plot(time_list,pz_error,label = r"$p_z$",ls = ":")
        ax[1].set_ylim(-0.1,0.1)
        ax[1].set_ylabel(r"Total Momentum Error")
        ax[1].legend(loc = 0)
        ax[2].plot(time_list,n_p_list)
        ax[2].set_ylabel("Number of Particles")
        ax[2].set_xlabel("Cycle number")
        fig.savefig(args["output"])
        plt.close(fig)
        if (mass_error_pass and px_error_pass and py_error_pass and pz_error_pass):
            print("All quantities conserved. Test passed")
            sys.exit(0)
        else:
            print("Not all quantities conserved. Test failed.")
            sys.exit(1)

    
