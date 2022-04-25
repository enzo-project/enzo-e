# ics.py generates a text file (called particles.dat by default)
# with particle data which can be used to initialise
# the merge stars test problem in Enzo-E. More specifically, it creates star
# particles with constant masses, with random positions uniformly distributed
# in a sphere, and initial velocities all of the same magnitude with directions
# directed towards the centre of the sphere, plus an additional (optional)
# constant drift velocity. For more details run "python ics.py -h".

import numpy as np
import argparse as ap

parser = ap.ArgumentParser(
    description="""
    Generates ICs for a sphere of star particles with random positions, on radial
    infalling tranjectories. 
    """
    )

parser.add_argument(
    "-r",
    "--radius",
    help="""
    The radius of the sphere in cm.
    """,
    required=True,
    type=float,
)


parser.add_argument(
    "-m",
    "--mass",
    help="The total mass of the sphere in grams.",
    required=True,
    type=float,
)

parser.add_argument(
    "-c",
    "--centre",
    nargs = 3,
    help="""
    The x,y,z coordinates of the centre of the sphere (in cm). Default = [0,0,0]
    """,
    required=False,
    type=float,
    default = [0.,0.,0.],
)

parser.add_argument(
    "-d",
    "--drift_velocity",
    nargs = 3,
    help="""
    The x,y,z components of the drift velocity of the sphere (in cm/s).
    Default = [0,0,0]
    """,
    required=False,
    type=float,
    default = [0.,0.,0.],
)

parser.add_argument(
    "-i",
    "--infall_speed",
    help="""
    The infall speed of the star particles in cm/s
    """,
    required=True,
    type=float,
)

parser.add_argument(
    "-f",
    "--filename",
    help="""
    The name of the file that is generated. Default = particles.dat
    """,
    required=False,
    default = "particles.dat",
    type=str,
)

parser.add_argument(
    "-s",
    "--seed",
    help="""
    The seed used to generate the particle positions. Default = 456
    """,
    required=False,
    default = 456,
    type=int,
)

parser.add_argument(
    "-n",
    "--n_particles",
    help="""
    The number of particles.
    """,
    required=True,
    type=int,
)

if __name__ == "__main__":
    args = vars(parser.parse_args())
    
    ## set seed
    np.random.seed(args["seed"])

    ## generate random numbers
    rand = np.random.rand(args["n_particles"],3)

    # (r,theta,phi polar coordinates)
    # Use first column to get r
    r = args["radius"] * np.cbrt(rand[:,0])
    # Use second column to get theta
    theta = np.arccos(1 - 2*rand[:,1])
    # Use third column to get phi
    phi = 2*np.pi*rand[:,2]

    # convert to cartesian coordinates
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)

    # get radial velocities
    vx = -1 * args["infall_speed"] * x / r
    vy = -1 * args["infall_speed"] * y / r
    vz = -1 * args["infall_speed"] * z / r

    # shift positions and velocities
    x += args["centre"][0]
    y += args["centre"][1]
    z += args["centre"][2]
    vx += args["drift_velocity"][0]
    vy += args["drift_velocity"][1]
    vz += args["drift_velocity"][2]

    # get masses
    mass = np.full(args["n_particles"],args["mass"]/args["n_particles"])

    # make the file
    output = np.column_stack((mass,x,y,z,vx,vy,vz))
    np.savetxt(args["filename"],output)
    
