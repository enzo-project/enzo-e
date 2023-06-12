from collections import namedtuple
import os
import os.path
import sys

import pytest

from answer_testing import set_cached_opts
from test_utils.enzoe_driver import EnzoEDriver

_root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))

# each entry in the following list specifies a command-line options that
# fall-back to an environment variable.
_ConfigParam = namedtuple(
    "_ConfigParam", ["flag", "env_var", "default", "help",
                     "other_argparse_kwargs"])
_CONFIG_OPTIONS = [
    _ConfigParam(
        flag = "--charm", env_var = "CHARM_PATH",
        default = None,
        help = ("path to charmrun binary that is used to execute enzo-e or "
                  "to the directory holding that binary."),
        other_argparse_kwargs = dict(action = "store")
    ),
    _ConfigParam(
        flag = "--local-dir", env_var = "TEST_RESULTS_DIR",
        default = "~/enzoe_test_results",
        help = "Path to directory where answers are/will be stored.",
        other_argparse_kwargs = dict(action = "store")
    ),
    _ConfigParam(
        flag = "--grackle-input-data-dir", env_var = "GRACKLE_INPUT_DATA_DIR",
        default = "~/enzoe_test_results",
        help = "Path to directory of grackle input data files.",
        other_argparse_kwargs = dict(action = "store")
    )
]


# this hook is used to add more command line flags to the pytest launcher
def pytest_addoption(parser):

    # add ordinary flags (NOT backed by an environment variable)

    # in the past, this could be set by the ENZO_PATH environment variable,
    # but this wasn't used by CI
    parser.addoption("--enzoe", action = 'store',
                     default = os.path.join(_root_dir, "build/bin/enzo-e"),
                     help = ("path to the enzo-e binary that will be tested. "
                             "Default is <ROOT_DIR>/build/bin/enzo-e"))


    # introduce flags that fall back to values specified in an env variable
    # -> the main reason these flags exist is for backwards compatability
    def _argparse_kwargs(other_kw, help, env_var):
        _ERR_PREFIX = ()
        if other_kw.get('action', 'store') in ['store_true', 'store_false']:
            raise RuntimeError(
                "for a flag backed by an environment variable, the argparse "
                "kwargs can't set 'action' to 'store_true' or 'store_false'. "
                "These flags must default to None")

        kwargs = other_kw.copy() # don't mutate the original
        kwargs["help"] = (f"{help} If not specified, the program tries to "
                          f"use the value set by {env_var}")
        kwargs["default"] = None
        return kwargs

    for param in _CONFIG_OPTIONS:
        parser.addoption(param.flag,
                         **_argparse_kwargs(param.other_argparse_kwargs,
                                            param.help, param.env_var))

    parser.addoption(
        "--answer-store",
        **_argparse_kwargs({'action' : "store_const", 'const' : True},
                           help="Indicates whether to generate test results.",
                           env_var = "GENERATE_TEST_RESULTS"))


def _to_abs_path(path):
    return os.path.abspath(os.path.expanduser(path))


# this hook gets executed after pytest finishes parsing all of the config opts
def pytest_configure(config):
    if "--help" in sys.argv[1:]:
        return

    # first, collect values from the command line
    vals = {'enzoe' : config.getoption('--enzoe')}

    for par in _CONFIG_OPTIONS: 
        arg_name = par.flag.lstrip('-').replace('-', '_')
        val = config.getoption(arg_name)
        if (val is None) and (os.environ.get(par.env_var, None) is None):
            if par.default is None:
                raise RuntimeError(
                    f"{par.arg_flag} or {par.env_var} must be set")
            else:
                val = par.default
        elif val is None:
            val = os.environ.get(par.env_var)
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
        uses_double_prec = enzoe_driver.query_uses_double_precision(),
        generate_results = vals['answer_store'],
        test_results_dir = _to_abs_path(vals['local_dir']),
        grackle_input_data_dir = _to_abs_path(vals['grackle_input_data_dir'])
    )
