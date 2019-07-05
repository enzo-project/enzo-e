import os

import numpy as np
import matplotlib.pyplot as plt
import yt

# the designed plots are provided to allow for easy comparison with figures
# 22 and 23 of Stone et al 2008 (note that their plots are for a higher order
# method)

def prep_cur_dir():
    cwd = os.getcwd()
    print(cwd[-22:])
    if cwd[-22:] == "input/vlct/orszag-tang":
        os.chdir("../../../")

def plot_panel(grid, field, title, ax, nlevels, idx = None, kwargs = {}):
    # field should either be a callable used to compute the data or is a
    # field name
    x = grid['x']
    y = grid['y']
    if callable(field):
        z = field(grid)
    else:
        z = grid[field]
    if idx is not None:
        x = x[idx]
        y = y[idx]
        z = z[idx]

    min_val = z.min()
    max_val = z.max()
    levels = np.linspace(min_val,max_val,num=nlevels+2)

    im = ax.contour(x,y,z, levels = levels[1:-1],**kwargs)
    ax.set_title(title)
    return im

def plot_ortho_ray(ax, ray, axis_name, field,fmt=None,kwargs = {}):
    data = ray[field]
    pos = ray[axis_name]

    idx = np.argsort(pos)

    data = data[idx]
    pos = pos[idx]

    if fmt is None:
        out = ax.plot(pos,data,**kwargs)
    else:
        out = ax.plot(pos,data,fmt,**kwargs)
    return out

if __name__ == '__main__':
    prep_cur_dir()

    ds = yt.load('method_vlct-8-orszag-tangN128_0.5/'
                 'method_vlct-8-orszag-tangN128_0.5.block_list')
    print(ds.parameters['Mesh']['root_size'])

    grid = ds.covering_grid(level = 0, left_edge = [0.,0.,0.],
                            dims=ds.domain_dimensions)


    # has 30 evenly spaced contours between minimum and maximum in each panel
    fig,ax_array = plt.subplots(2,2,sharex=True,sharey=True, figsize = (6,6))

    kwargs = {'linewidths' : 1}

    plot_panel(grid, "density", "density", ax_array[0,0], nlevels=30,
               idx = (slice(None),slice(None),0),
               kwargs = kwargs)

    plot_panel(grid, "pressure", "pressure", ax_array[0,1], nlevels=30,
               idx = (slice(None),slice(None),0),
               kwargs = kwargs)

    ax_array[0,1].axhline(0.3125,color = 'k', linestyle = '--', linewidth = 1)
    ax_array[0,1].axhline(0.427,color = 'k', linestyle = '--', linewidth = 1)

    plot_panel(grid, "magnetic_pressure", "magnetic pressure", ax_array[1,0],
               nlevels=30, idx = (slice(None),slice(None),0),
               kwargs = kwargs)

    calc_specific_ke = lambda g : g["kinetic_energy"]/g["density"]
    plot_panel(grid, calc_specific_ke, "specific kinetic energy",
               ax_array[1,1], nlevels=30, idx = (slice(None),slice(None),0),
               kwargs = kwargs)

    for ax in ax_array.flatten():
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.tick_params(axis="y",direction="in")
        ax.tick_params(axis="x",direction="in")
        ax.set_aspect('equal')

    fig.tight_layout()
    plt.savefig("input/vlct/orszag-tang/0.5_contours.pdf")
    plt.close('all')
    

    fig,ax_array = plt.subplots(2,1,figsize=(6,4),sharex=True)
    kwargs = {'mfc' : 'none'}

    oray = ds.ortho_ray("x", (0.3125, 0.0))
    plot_ortho_ray(ax_array[0], oray, 'x', 'pressure', fmt='ks',
                   kwargs = kwargs)

    oray = ds.ortho_ray("x", (0.427, 0.0))
    plot_ortho_ray(ax_array[1], oray, 'x', 'pressure', fmt='ks',
                   kwargs = kwargs)
    ax_array[1].set_xlim(0.,1.)
    ax_array[0].set_ylabel('P')
    ax_array[1].set_ylabel('P')

    ax_array[0].set_ylim(0.0,0.3)
    ax_array[1].set_ylim(0.0,0.4)
    ax_array[1].set_xlabel('x')

    fig.tight_layout()

    plt.savefig("input/vlct/orszag-tang/pressure_slice.pdf")

