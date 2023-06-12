from collections import namedtuple
from dataclasses import asdict
import os
import os.path
import sys

import pytest

from answer_testing import set_cached_opts, cached_opts
from test_utils.enzoe_driver import EnzoEDriver
from test_utils.parse_cmake_cache import parse_cmake_cache

# all cmd flags that fall-back to an environment variable are listed below
_ConfigParam = namedtuple("_ConfigParam", ["env_var", "default", "help",
                                           "other_argparse_kwargs",
                                           "coerce_env_val"])
_CONFIG_OPTIONS = {
    '--charm' : _ConfigParam(
        env_var = "CHARM_PATH",
        default = None,
        help = ("path to charmrun binary that is used to execute enzo-e or "
                "to the directory holding that binary."),
        other_argparse_kwargs = dict(action = "store"),
        coerce_env_val = lambda s: s
    ),
    "--local-dir" : _ConfigParam(
        env_var = "TEST_RESULTS_DIR",
        default = "~/enzoe_test_results",
        help = "Path to directory where answers are/will be stored.",
        other_argparse_kwargs = dict(action = "store"),
        coerce_env_val = lambda s: s
    ),
    "--grackle-input-data-dir" : _ConfigParam(
        env_var = "GRACKLE_INPUT_DATA_DIR",
        default = None,
        help = "Path to directory of grackle input data files.",
        other_argparse_kwargs = dict(action = "store"),
        coerce_env_val = lambda s: s
    ),
    # it's important that the following is store_const instead of store_true
    # (so that we can force argparse value default to None)
    "--answer-store" : _ConfigParam(
        env_var = "GENERATE_TEST_RESULTS",
        default = False,
        help = "Indicates whether to generate test results.",
        other_argparse_kwargs = dict(action = "store_const", const = True),
        coerce_env_val = lambda s: s.lower() == "true"
    )
}


_root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))


# this hook is used to add more command line flags to the pytest launcher
def pytest_addoption(parser):

    # add ordinary flags (NOT backed by an environment variable)

    # in the future, it would probably be more ergonomic to assume that the
    # current working directory is the enzo-e build directory when this flag is
    # not specified (and we may want to try parsing CMakeCache.txt)
    parser.addoption(
        "--build-dir", action = 'store', default = None,
        help = ("points to the build-directory where the target enzo-e binary "
                "was built (that binary has the path: BUILD_DIR/bin/enzo-e). "
                "The path to the charmrun launcher will be inferred from the "
                "BUILD_DIR/CMakeCache.txt file, but can be overwritten by the "
                "--charm flag or the CHARM_PATH environment variable. Be "
                "aware that problems can theoretically arise if you do "
                "something to alter the build system's state after a build "
                "without rebuilding enzo-e (CMakeCache.txt may no longer be "
                "valid for the binary). When this flag isn't specified, the "
                "test-driver assumes that the enzo-e binary is located at "
                f"{_root_dir}/build/bin/enzo-e, but does not try to infer "
                "charmrun's location from CMakeCache.txt.")
    )


    # introduce flags that fall back to values specified in an env variable
    # (these flags primarily exist for backwards compatability)
    for flag, par in _CONFIG_OPTIONS.items():
        kwargs = par.other_argparse_kwargs.copy() # don't mutate the original
        kwargs["help"] = (f"{par.help} If not specified, the program tries to "
                          f"use the value set by {par.env_var}")
        kwargs["default"] = None
        if kwargs.get('action', 'store') in ['store_true', 'store_false']:
            raise RuntimeError(
                "for a flag backed by an environment variable, the argparse "
                "kwargs can't set 'action' to 'store_true' or 'store_false'.")
        parser.addoption(flag, **kwargs)


def _to_abs_path(path):
    if path is None:
        return None
    return os.path.abspath(os.path.expanduser(path))


# this hook gets executed after pytest finishes parsing all of the config opts
def pytest_configure(config):
    # inspect/parse parameters and initialize the globally cached options

    # step 1: early exit if --help was passed
    if "--help" in sys.argv[1:]:
        return

    # step 2: fetch the parameters that are backed by environment variables
    vals = {}
    for flag, par in _CONFIG_OPTIONS.items(): 
        arg_name = flag.lstrip('-').replace('-', '_')
        if config.getoption(arg_name) is not None:
            vals[arg_name] = config.getoption(arg_name)
        elif os.environ.get(par.env_var, None) is not None:
            vals[arg_name] = par.coerce_env_val(os.environ.get(par.env_var))
        else:
            vals[arg_name] = par.default
    
    # step 3: initialize the enzo-e driver
    build_dir = _to_abs_path(config.getoption("build_dir"))

    # step 3a: determine the path to the enzo-e binary
    if build_dir is not None:
        enzoe_path = os.path.join(build_dir, 'bin/enzo-e')
        if not os.path.isfile(enzoe_path):
            raise RuntimeError(f"enzo-e isn't at {enzoe_path}.")
    else:
        enzoe_path = os.path.join(_root_dir, 'build/bin/enzo-e')
        if not os.path.isfile(enzoe_path):
            raise RuntimeError(f"enzo-e isn't at {enzoe_path}. If the it is "
                               "located elsewherer use --build-dir.")

    # step 3b: determine the path to charmrun
    if (vals['charm'] is None) and (build_dir is None):
        raise RuntimeError("Please specify --charm or set CHARM_PATH when "
                           "--build-dir isn't specified")
    elif vals['charm'] is None:
        # this approach may be brittle (which is why we give --charm priority)
        # - if it is brittle, we could have cmake directly spit out a JSON file 
        #   with the desired info as part of the build 
        _cache = parse_cmake_cache(os.path.join(build_dir, "CMakeCache.txt"))
        if ("CHARM_RUN" not in _cache) or (_cache["CHARM_RUN"] == ''):
            raise RuntimeError("did the buildsystem change? CHARM_RUN was not "
                               "cached or it was assigned an empty string")
        charmrun_path = _cache["CHARM_RUN"]
    else:
        charmrun_path = _to_abs_path(vals['charm'])
        if os.path.basename(charmrun_path) != "charmrun":
            charmrun_path = os.path.join(charmrun_path, 'charmrun')

    # Step 3c: initialize EnzoEDriver
    enzoe_driver = EnzoEDriver(enzoe_path = enzoe_path,
                               charmrun_path = charmrun_path)

    # Step 4: initialize the cached options
    set_cached_opts(
        enzoe_driver = enzoe_driver,
        uses_double_prec = enzoe_driver.query_uses_double_precision(),
        generate_results = vals['answer_store'],
        test_results_dir = _to_abs_path(vals['local_dir']),
        grackle_input_data_dir = _to_abs_path(vals['grackle_input_data_dir'])
    )


# hook for printing information at the top of the test report
def pytest_report_header(config):
    opts = cached_opts()

    def format_kv_pairs(prefix, pairs):
        blank_prefix = ' ' * len(prefix)
        for i, (k,v) in enumerate(pairs):
            if i == 0:
                yield f'{prefix}{k} = {v}'
            else:
                yield f'{blank_prefix}{k} = {v}'

    part1 = list(format_kv_pairs("enzoe_driver:     ",
                                 asdict(opts.enzoe_driver).items()))
    part2 = list(format_kv_pairs(
        "other properties: ",
        filter(lambda pair: pair[0] != "enzoe_driver", asdict(opts).items())
    ))

    return part1 + part2
