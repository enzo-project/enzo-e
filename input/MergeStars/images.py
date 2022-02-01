import yt
import matplotlib.pyplot as plt
import argparse as ap
import matplotlib
matplotlib.use("Agg")

parser = ap.ArgumentParser(
    description="""
    Reads a data directory containing a series of snapshots of the MergeStarsTest
    problem and plots the x and y positions of the star particles, making a separate
    image for each snapshot
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
args = vars(parser.parse_args())
input_prefix = args["input_prefix"]
output_prefix = args["output_prefix"]
yt.enable_parallelism()
ds_pattern = f"{input_prefix}_????/{input_prefix}_????.block_list"
ts =  yt.DatasetSeries(ds_pattern)

storage = {}
for ds in ts.piter():

    box = ds.box(left_edge = -ds.domain_width/2.0,right_edge = ds.domain_width/2.0)
    x = box["star","x"]
    y = box["star","y"]
    z = box["star","z"]
    fig,ax = plt.subplots()
    ax.plot(x,y,marker = "x",linestyle = "None")
    ax.set_xlim(-ds.domain_width[0]/2.0,ds.domain_width[0]/2.0)
    ax.set_ylim(-ds.domain_width[0]/2.0,ds.domain_width[0]/2.0)
    ax.set_xlabel("x (cm)")
    ax.set_ylabel("y (cm)")
    current_cycle = ds.parameters["current_cycle"]
    filename = f"{output_prefix}_{current_cycle}.png"
    fig.savefig(filename)
    plt.close(fig)

    
