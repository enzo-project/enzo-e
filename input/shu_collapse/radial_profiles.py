#!/bin/python

# This file defines a function called make_radial_profiles, which
# takes as input the directories generated when Enzo-E runs the test
# problem, which contain snapshots of the particle data and field data
# at regular time intervals.
# It then uses yt make plots of spherically averaged quantities (density, radial velocity,
# and mass flux) as a function of radius, at times specified by `profile_times_years`,
# with the coordinates of the center of the sphere specified by `center`. The plots are saved
# in a file called `radial_profiles.png`.
# These plots can be compared to (for example) Figure 10 of Federrath et al 2010, ApJ, 713, 269.
# make_radial_profiles can be imported by another script (as is done in
# run_shu_collapse.py), or it can be executed by running this file
# as a script, with command line arguments being passed to make_radial_profiles. 
# For more information, run "python radial_profiles.py -h".

import yt
from yt.fields.api import ValidateParameter
import numpy as np
import matplotlib.pyplot as plt
import glob
import unyt
import argparse as ap
import sys

def make_radial_profiles(center,profile_times_years,input_prefix):

    def _my_radial_velocity(field,data):
        xv = data["gas","velocity_x"]
        yv = data["gas","velocity_y"]
        zv = data["gas","velocity_z"]
        x_hat = data["enzoe","x"] - center[0]
        y_hat = data["enzoe","y"] - center[1]
        z_hat = data["enzoe","z"] - center[2]
        r = np.sqrt(x_hat*x_hat+y_hat*y_hat+z_hat*z_hat)
        x_hat /= r
        y_hat /= r
        z_hat /= r
        return xv*x_hat + yv*y_hat + zv*z_hat

    def _mass_flux(field,data):

        xv = data["gas","velocity_x"]
        yv = data["gas","velocity_y"]
        zv = data["gas","velocity_z"]
        x_hat = data["enzoe","x"] - center[0]
        y_hat = data["enzoe","y"] - center[1]
        z_hat = data["enzoe","z"] - center[2]
        r = np.sqrt(x_hat*x_hat+y_hat*y_hat+z_hat*z_hat)
        x_hat /= r
        y_hat /= r
        z_hat /= r
        rad_vel = xv*x_hat + yv*y_hat + zv*z_hat

        density = data["gas","density"]

        return density*rad_vel

    ds_pattern = f"{input_prefix}_????/{input_prefix}_????.block_list"
    time_series = yt.load(ds_pattern)
    time_list = [ds.current_time.to("yr") for ds in time_series]
    ind_list = [find_snapshot_index(t,time_list) for t in profile_times_years]

    yt.add_field(("gas","my_radial_velocity"),
             sampling_type = "cell",
             function=_my_radial_velocity,
             units="cm/s",
             take_log=False,
             validators=[ValidateParameter(['center'])])
    yt.add_field(("gas","mass_flux"),
             sampling_type = "cell",
             function=_mass_flux,
             units="g/cm**2/s",
             take_log=False,
             validators=[ValidateParameter(['center'])])

    profiles = [get_profiles(time_series,time_list,ind,center)
                for ind in ind_list if ind >= 0]

    fig,ax = plt.subplots(nrows = 3,sharex = True)
    fig.set_figheight(12)
    fig.set_figwidth(6)

    for profile in profiles:
        ax[0].plot(profile["radius"],profile["density"])
        ax[0].set_ylim(1.0e-19,1.0e-14)
        ax[0].set_ylabel(r"$\rho\,(g\,cm^{-3})$")
        ax[0].set_yscale('log')
        ax[0].axvline(x = 4.0*profile["cell_width"],ls = "--",color = "k",lw = 0.5)
        ax[1].plot(profile["radius"],profile["radial_velocity"])
        ax[1].set_ylim(-5,1)
        ax[1].set_ylabel(r"$v_r \, (km\,s^{-1}$")
        ax[1].axvline(x = 4.0*profile["cell_width"],ls = "--",color = "k",lw = 0.5)
        ax[2].plot(profile["radius"],profile["mass_flux"],
                   label = r"$%1.1f \times 10^3\,\mathrm{yrs}$" %(profile["time"]/1.0e3))
        ax[2].set_yscale('log')
        ax[2].set_ylim(1.0e-6,1.0e-2)
        ax[2].set_ylabel(r"$\dot{M} \, (M_{\odot}/\mathrm{yr})$")
        ax[2].set_xlabel(r"$r\,cm$")
        ax[2].set_xlim(1.0e15,1.0e17)
        ax[2].set_xscale('log')
        ax[2].axvline(x = 4.0*profile["cell_width"],ls = "--",color = "k",lw = 0.5)
        ax[2].axhline(y = 1.5e-4,ls = "--",color = "k")
    ax[2].legend(bbox_to_anchor = (0.0,-0.2),loc = "upper left",ncol = 3)

    fig.savefig("radial_profiles.png")

def find_snapshot_index(target_time,time_list):

    # finds closest time in time_list to the target time by using binary search
    # assumes time_list in sorted
    if target_time > time_list[-1]:
        print(f"Warning: target_time = {target_time}, largest value in time_list is "
              f"{time_list[-1]:3g}")
        return -1
    elif target_time < time_list[0]:
        print(time_list[0])
        print(f"Warning: target_time = {target_time}, smallest value in time_list is "
              f"{time_list[0]}:3g")
        return -1

    else:
        lower_bound = 0
        upper_bound = len(time_list) - 1
        while (upper_bound - lower_bound > 1):
            midpoint = int((lower_bound+upper_bound)/2)
            if (target_time >= time_list[midpoint]):
                lower_bound = midpoint
            else:
                upper_bound = midpoint

        if target_time < 0.5*(time_list[lower_bound] + time_list[upper_bound]):
            return lower_bound
        else:
            return upper_bound

def get_profiles(time_series,time_list,ind,center):

    data = {}
    ds = time_series[ind]
    sphere = ds.sphere(center,ds.domain_width[0]/2.0)
    profiles = yt.create_profile(sphere,"radius",["density","my_radial_velocity","mass_flux"],
                                 weight_field= "cell_volume",logs = {"radius":False},
                                   units = {'radius': 'cm',
                                            'density':"g/cm**3",
                                            'my_radial_velocity':'km/s',
                                            "mass_flux": "g/cm**2/s"})
    data["time"] = time_list[ind]
    data["radius"] = profiles.x.value
    data["density"] = profiles["density"].value
    data["radial_velocity"] = profiles["my_radial_velocity"].value
    mass_flux  = -4.0 * np.pi * profiles.x**2 * profiles["mass_flux"]
    data["mass_flux"] = mass_flux/(unyt.M_sun/unyt.year)
    data["cell_width"] = ds.domain_width[0]/ds.domain_dimensions[0]

    np.clip(data["mass_flux"],a_min = 0.0,a_max = None)

    return data

if __name__ == "__main__":

    parser = ap.ArgumentParser()
    parser.add_argument(
        "-c",
        "--center",
        nargs = 3,
        help="""
        The x,y,z coordinates of the center of the collapsing sphere (in code units).
        """,
        required=True,
        type=float,
        default = [0.,0.,0.])

    parser.add_argument(
        "-t",
        "--profile_times_years",
        nargs = 3,
        help="""
        The times (in years) at which radial profiles are computed.
            Default = [0.0,8.0e3,1.2e4,1.5e4,1.7e4].
        """,
        required=False,
        type=float,
        default = [0.0,8.0e3,1.2e4,1.5e4,1.7e4])


    parser.add_argument(
        "-i",
        "--input_prefix",
        help="The prefix for the directories containing snapshots. Default = Dir",
        required=False,
        default = "Dir",
        type=str,
    )

    args = parser.parse_args()

    make_radial_profiles(args.center,args.profile_times_years,args.input_prefix)
