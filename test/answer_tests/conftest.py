import os
import os.path

import pytest

from answer_testing import set_cached_opts
from test_utils.enzoe_driver import EnzoEDriver

# each entry in this list is a tuple with the following format:
# (name or flags -- to be passed to argparse.add_arguments,
#  kwargs to pass to argparse.add_arguments,
#  fallback-env-variable,
#  fallback_value -- when neither the flag nor env-var is provided)
_CONFIG_OPTIONS = [
    ('--enzoe',
     dict(action = "store", default = None,
          help = "path to the enzo-e binary that will be tested."),
     "ENZO_PATH", "build/bin/enzo-e"),

    ("--charm",
     dict(action="store", default = None,
          help = ("path to charmrun binary that is used to execute enzo-e or "
                  "to the directory holding that binary.")),
     "CHARM_PATH", None),

    ("--local-dir",
     dict(action="store", default = None,
          help = "Path to directory where answers are/will be stored."),
     "TEST_RESULTS_DIR", "~/enzoe_test_results"),

    ("--grackle-input-data-dir",
     dict(action="store", default = None,
          help = "Path to directory of grackle input data files."),
     "GRACKLE_INPUT_DATA_DIR", None),
]


# this hook is used to add more tests to the pytest launcher
def pytest_addoption(parser):

    def _fallback_help_clause(env_var):
        return ("If not specified, the program tries to use the value set by "
                f"the {env_var} environment variable (if specified)")

    for arg_flag, kwargs, env_var, _ in _CONFIG_OPTIONS:
        tmp = kwargs.copy()
        tmp["help"] += " " + _fallback_help_clause(env_var)
        parser.addoption(arg_flag, **tmp)

    parser.addoption(
        "--answer-store", action="store_const",
        default = None, const = True,
        help = ("Indicates whether to generate test results. " +
                _fallback_help_clause("GENERATE_TEST_RESULTS"))
    )


def _to_abs_path(path):
    return os.path.abspath(os.path.expanduser(path))


# this hook gets executed after pytest finishes parsing all of the config opts
def pytest_configure(config):

    # first, collect values from the command line
    vals = {}
    for arg_flag, _, env_var, fallback in _CONFIG_OPTIONS:
        arg_name = arg_flag.lstrip('-').replace('-', '_')
        val = config.getoption(arg_name)
        if (val is None) and (os.environ.get(env_var, None) is None):
            if fallback is None:
                raise RuntimeError(f"{arg_flag} or {env_var} must be set")
            else:
                val = fallback
        else:
            val = os.environ.get(env_var)
        vals[arg_name] = val

    if config.getoption('answer_store') is None:
        vals['answer_store'] = os.environ.get(
            "GENERATE_TEST_RESULTS", "false").lower() == "true"
    else:
        vals['answer_store'] = config.getoption('answer_store')

    # now, lets initialize the enzo-e driver
    charmrun_path = _to_abs_path(vals['charm'])
    if os.path.basename(charmrun_path) != "charmrun":
        charmrun_path = os.path.join(charmrun_path, 'charmrun')
    enzoe_driver = EnzoEDriver(enzoe_path = _to_abs_path(vals['enzoe']),
                               charmrun_path = charmrun_path)

    set_cached_opts(
        enzoe_driver = enzoe_driver,
        generate_results = vals['answer_store'],
        test_results_dir = _to_abs_path(vals['local_dir']),
        grackle_input_data_dir = _to_abs_path(vals['grackle_input_data_dir'])
    )
