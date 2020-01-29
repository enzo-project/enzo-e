import sys
import math

import numpy as np

# Determining the velocity:
# In the normal case, we use a domain with stretching from 0 to 1 with 256
# cells; this gives a dx = 1. / 256 . We also know that the interval of time
# that we allow the cloud to evolve over is dt = 0.25 . Additionally, we can
# compute the sound speed in the hottest gas of the initial conditions using
# cs = sqrt (gamma * p / rho); canonically, gamma=1.4, p = 1, and rho = 1.
#
# If we are targetting some Mhat then we can determine the closest velocity
# that leads to a translation over an integral number of cells over dt. The
# number of cells is given by:
#    Ncells = floor(cs * Mhat / (dx / dt) + 0.5)
# We could further require that Ncells be divisible by 4:
#    Ncells = floor(floor(cs * Mhat / (dx / dt) + 0.5)/4 + 0.5)
# Once we know Ncells, we can easily compute the additional velocity component
#    velocity = Ncells * dx / dt

# the interval over which the problem usually evolves (without the
# transformation)
_DEFAULT_XLIM = (0.,1.0)

# resolution of the cells
_ORIG_CELLS = 128

# time over which the system evolves
_DT = 0.25

# designates that x-axis must have a total number of cells divisable by 4:
_TOTCELLS_REQ_FACTOR = 4

# adiabatic index
_GAMMA = 1.4

if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise RuntimeError("A single command line argument is required to "
                           "denote the target Mach number")
    mach_number_arg = float(sys.argv[-1])

    pos_vel = (mach_number_arg > 0)
    mach_target = math.fabs(mach_number_arg)

    cs = math.sqrt(_GAMMA*1.0/1.0)

    default_interval = _DEFAULT_XLIM[1] - _DEFAULT_XLIM[0]
    dx = default_interval/float(_ORIG_CELLS)

    # compute number of extra cells - this needs to be adjusted to make sure
    # the total number of cells is divisible by _TOTCELLS_REQ_FACTOR
    n_extra_cells_guess = math.floor((cs * mach_target * _DT / dx) + 0.5)
    tot_cells_guess = _ORIG_CELLS + n_extra_cells_guess

    tot_cells = int(math.floor(tot_cells_guess/float(_TOTCELLS_REQ_FACTOR)
                               + 0.5) * _TOTCELLS_REQ_FACTOR)
    n_extra_cells = tot_cells - _ORIG_CELLS

    velocity = n_extra_cells * dx / _DT
    if not pos_vel:
        velocity*=-1.
    final_mach = velocity / cs

    msg = "target mach: {}, final mach: {}, final velocity: {}"

    print(msg.format(repr(mach_target),repr(final_mach),
                     repr(velocity)))

    # now, provide the new domain
    if pos_vel:
        new_lower = _DEFAULT_XLIM[0]
        new_upper = tot_cells*dx + new_lower

        extract_cells = [n_extra_cells, tot_cells]
        extract_range = [n_extra_cells*dx + new_lower,new_upper]
        
    else:
        new_upper = _DEFAULT_XLIM[1]
        new_lower = new_upper - tot_cells*dx
        extract_cells = [0, _ORIG_CELLS]
        extract_range = [new_lower,new_upper-n_extra_cells*dx]

    msg = "res = {!r}/{!r}, tot_cells: {!r}, new_interval: [{!r}, {!r}]"
    print(msg.format(default_interval, _ORIG_CELLS, tot_cells,
                     new_lower, new_upper))

    # finally print out the interval from which the solution should be extracted
    msg = ("check results in cells with (centered) positions between: "
           "[{!r},{!r}]")
    print(msg.format(*extract_range))
    msg = "these active zone indices, start at {!r} and stop at {!r}"
    print(msg.format(*extract_cells))
