from math import ceil
import os
import os.path
import yt
import libconf # for help reading the parameter file (dependency of yt)

from answer_testing import EnzoETest

_base_file = os.path.basename(__file__)

def _sched_var_type(name):
    d = {"cycle" : int, "time" : float,
         "seconds" : float, "minutes" : float, "hours" : float}
    try:
         return d[name]
    except KeyError:
        raise ValueError(
            f"{name!r} is not a known scheduling variable. Those include "
            f"{list(name_type_pairs.keys())!r}."
        )

def _all_schedule_var_vals(sched_params, var_name, start_val, last_val):
    """
    Returns a sequence holding all values of a schedule variable in the 
    inclusive interval [start_val, last_val]
    """

    if sched_params.get('var', None) != var_name:
        raise ValueError("This function is expecting for there to be a 'var' "
                         f"parameter in sched_params with value {var_name!r}")

    var_type = _sched_var_type(var_name)
    def coerce_arg_(arg_name, arg):
        out = var_type(arg)
        if out != arg:
            raise ValueError(
                f"the argument {arg_name!r} has a value {arg} that can't be "
                f"losslessly coerced to the {var_type.__name__} dtype "
                f"associated with the {var_name!r} scheduling variable."
            )
        return out

    start_val = coerce_arg_("start_val", start_val)
    last_val = coerce_arg_("last_val", last_val)
    assert last_val > start_val

    if ('list' in sched_params) and (len(sched_params) == 2):
        return sorted([e for e in sched_params['list']
                       if ((e >= start_val) and (e <= last_val))])
    elif all(k in ('var', 'start', 'stop', 'step') for k in sched_params):
        # parse parameters
        sched_start = sched_params.get('start', 0)
        sched_step = sched_params.get('step', 1)
        if (int(sched_start)!=sched_start) or (int(sched_step)!=sched_step):
            # we can run into problems with exact floating point values if we
            # don't do this...
            raise RuntimeError("This function only handles schedule:start and "
                               "schedule:step values that are integers")

        # the schedule:stop parameter is inclusive, so fall back to last_val
        # (also inclusive) when schedule:stop isn't specified
        sched_stop = sched_params.get('stop', last_val)

        # compute vals
        if sched_stop < start_val:
            return []
        cur_val = sched_start
        # if applicable, advance cur_val until, cur_val >= start_val
        while cur_val <= start_val:
            cur_val += sched_step
        out = []
        while cur_val <= min(last_val, sched_stop):
            out.append(cur_val)
            cur_val += sched_step
        return out
    else:
        raise ValueError(f"{sched_params!r} doesn't describe a valid schedule")

def scheduled_output_dirs(parameters, var_name, start_val, last_val):
    """
    Constructs a generator that iterates over the expected output directories 
    produced for a given configuration of a given output-method that uses a
    schedule parameterized by the scheduling-variable, var_name.

    Specifically, this iterates over directories when the scheduling-variable
    has values in the inclusive interval [start_val, last_val]. This yields a
    2-tuple where the first element gives the directory name and the second
    gives the associated value of the scheduling-variable.

    Notes
    -----
    When the scheduling-variable is associated with a floating point variable
    (e.g. ``var_name == "time"``), one could imagine that this generator might
    slightly mis-predict the exact values of the scheduling-variable when the
    output-directories are produced. For that reason, this generator requires
    that the schedule:start and schedule:stop parameters are integers, when
    they are specified.

    This is a fairly generic operation that could be helpful in a lot of
    scenarios.
    """

    # construct a function to format the path_names parameter for a given value
    # of the scheduling parameter
    if ('path_name' not in parameters):
        raise ValueError("output-method is missing the path_name parameter")
    elif not isinstance(parameters['path_name'], (list,tuple,str)):
        raise ValueError("path_name parameter has unexpected type")
    elif isinstance(parameters['path_name'], str):
        def format_path(val): return parameters['path_name']
    elif len(parameters['path_name']) == 1:
        def format_path(val): return parameters['path_name'][0]
    elif any(fmt_v != var_name for fmt_v in parameters['path_name'][1:]):
        # we could relax this requirement if we knew the exact mapping between
        # the various scheduling variables (i.e. from parsing the log)
        raise ValueError("{var_name!r} is expected to be the only variable "
                         "used for formatting the path_name parameter")
    else:
        def format_path(val):
            n_vars = len(parameters['path_name']) - 1
            fmt_arg = tuple(val for i in range(n_vars))
            return parameters['path_name'][0] % fmt_arg

    if (('schedule' not in parameters) or
        (parameters['schedule']['var'] != var_name)):
        raise ValueError("output-method must have an associated schedule that "
                         f"is parameterized by the {var_name!r} variable.")
    vals = _all_schedule_var_vals(parameters['schedule'], var_name,
                                  start_val, last_val)

    for val in vals:
        yield format_path(val), val


class TestMethodScheduleTime(EnzoETest):
    """
    The purpose of this test is to ensure that the Enzo-E executes methods at 
    the exact simulation times at when they are scheduled.

    In short, we setup a minimal simulation and use the "null" Method to 
    enforce an upper-bound on the timestep during each cycle. We also setup 2
    occurences of the output method with schedules based on simulation-times, 
    that include execution times that shouldn't coincide occur when just
    considering the maximum timestep. This python code checks when the outputs
    are found.

    Unlike many other tests - this test does not require comparison against 
    previous versions of the code. There is an unambiguous "right answer" that
    shouldn't change...
    """

    parameter_file = "Schedule/Method_Scheduling_time.in"
    max_runtime = 5
    ncpus = 1

    def test_method_schedule_time(self):
        # parse the parameters (after they were translated to libconfig format)
        with open('parameters.libconfig', 'r') as f:
            config = libconf.load(f)

        start_t = 0.0
        stop_t = config["Stopping"]["time"]

        itr = scheduled_output_dirs(config["Method"]["output_step"],
                                    "time", start_t, stop_t)
        for dirname, t in itr:
            if t != stop_t:
                assert os.path.isdir(dirname)

        itr = scheduled_output_dirs(config["Method"]["output_list"],
                                    "time", start_t, stop_t)
        for dirname, t in itr:
            if t != stop_t:
                assert os.path.isdir(dirname)
