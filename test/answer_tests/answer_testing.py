import copy
import numpy as np
import os
import pytest
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import yt

from numpy.testing import assert_array_equal
from unittest import TestCase
from yt.funcs import ensure_dir
from yt.testing import assert_rel_equal

_base_file = os.path.basename(__file__)

# If GENERATE_TEST_RESULTS="true", just generate test results.
generate_results = os.environ.get("GENERATE_TEST_RESULTS", "false").lower() == "true"
yt.mylog.info(f"{_base_file}: generate_results = {generate_results}")

_results_dir = os.environ.get("TEST_RESULTS_DIR", "~/enzoe_test_results")
test_results_dir = os.path.abspath(os.path.expanduser(_results_dir))
yt.mylog.info(f"{_base_file}: test_results_dir = {test_results_dir}")
if generate_results:
    ensure_dir(test_results_dir)
else:
    if not os.path.exists(test_results_dir):
        raise RuntimeError(
            f"Test results directory not found: {test_results_dir}.")

# Set the path to charmrun
_charm_path = os.environ.get("CHARM_PATH", "")
if not _charm_path:
    raise RuntimeError(
        f"Specify path to charm with CHARM_PATH environment variable.")
charmrun_path = os.path.join(_charm_path, "charmrun")
yt.mylog.info(f"{_base_file}: charmrun_path = {charmrun_path}")
if not os.path.exists(charmrun_path):
    raise RuntimeError(
        f"No charmrun executable found in {_charm_path}.")

src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))

# Set the path to the enzo-e binary
_enzo_path = os.environ.get("ENZO_PATH", "")
if _enzo_path:
    enzo_path = os.path.abspath(_enzo_path)
else:
    enzo_path = os.path.join(src_path, "build/bin/enzo-e")
yt.mylog.info(f"{_base_file}: enzo_path = {enzo_path}")
if not os.path.exists(enzo_path):
    raise RuntimeError(
        f"No enzo-e executable found in {enzo_path}.")

input_dir = "input"

# check for data path to grackle
_grackle_input_data_dir = os.environ.get("GRACKLE_INPUT_DATA_DIR", None)
if ((_grackle_input_data_dir is not None) and
    (not os.path.exists(_grackle_input_data_dir))):
    raise RuntimeError("GRACKLE_INPUT_DATA_DIR points to a non-existent "
                       f"directory: {_grackle_input_data_dir}")


_grackle_tagged_tests = set()

def uses_grackle(cls):
    """
    Decorator that annotates that a test class uses grackle

    In detail, this annotates and sets up the appropriate skipif marker
    and updates a global registry of tests using grackle.
    """
    _grackle_tagged_tests.add(cls.__name__)

    wrapper_factory = pytest.mark.skipif(
        _grackle_input_data_dir is None,
        reason = "GRACKLE_INPUT_DATA_DIR is not defined"
    )
    return wrapper_factory(cls)

class EnzoETest(TestCase):
    parameter_file = None
    max_runtime = np.inf
    ncpus = None

    def setup_symlinks(self):
        ipath = os.path.join(src_path, input_dir)
        spath = os.path.join(self.tmpdir, input_dir)
        os.symlink(ipath, spath)

        if self.__class__.__name__ in _grackle_tagged_tests:
            # make symlinks to each grackle input file
            with os.scandir(_grackle_input_data_dir) as it:
                for entry in it:
                    if not entry.is_file():
                        continue
                    os.symlink(entry.path,os.path.join(self.tmpdir, entry.name))

    def run_simulation(self):
        pfile = os.path.join(input_dir, self.parameter_file)
        command = f"{charmrun_path} ++local +p{self.ncpus} {enzo_path} {pfile}"
        proc = subprocess.Popen(
            command, shell=True, close_fds=True,
            preexec_fn=os.setsid)

        stime = time.time()
        while proc.poll() is None:
            if time.time() - stime > self.max_runtime:
                os.killpg(proc.pid, signal.SIGUSR1)
                raise RuntimeError(
                    f"Simulation {self.__class__.__name__} exceeded max runtime of "
                    f"{self.max_runtime} seconds.")
            time.sleep(1)

        if proc.returncode != 0:
            raise RuntimeError(
                f"Simulation {self.__class__.__name__} exited with nonzero return "
                f"code {proc.returncode}.")

    def setUp(self):
        self.curdir = os.getcwd()
        self.tmpdir = tempfile.mkdtemp()

        self.setup_symlinks()
        os.chdir(self.tmpdir)
        self.run_simulation()

    def tearDown(self):
        os.chdir(self.curdir)
        shutil.rmtree(self.tmpdir)

def ytdataset_test(compare_func, **kwargs):
    """
    ytdataset test decorator.

    Put this decorator above testing functions that return a
    dictionary of values.
    """

    def real_answer_test(func):
        def wrapper(*args):
            # name the file after the function
            filename = "%s.h5" % func.__name__
            result_filename = os.path.join(test_results_dir, filename)

            # check that answers exist
            if not generate_results:
                assert os.path.exists(result_filename), \
                  "Result file, %s, not found!" % result_filename

            data = func(*args)
            fn = yt.save_as_dataset({}, filename=filename, data=data)

            # if generating, move files to results dir
            if generate_results:
                shutil.move(filename, result_filename)
            # if comparing, run the comparison
            else:
                ytdataset_compare(
                    filename, result_filename,
                    compare_func=compare_func, **kwargs)
        return wrapper
    return real_answer_test

def assert_array_rel_equal(a1, a2, decimals=16, **kwargs):
    """
    Wraps assert_rel_equal with, but decimals is a keyword arg.

    Adds functionality to interpret `decimals = np.inf` as an indication that
    an exact match is required
    """
    if decimals == np.inf:
        diff = a1-a2 # if incompatible units are attached to a1 and a2, an
                     # exception will be raised
        if hasattr(diff, 'ndarray_view'):
            diff = diff.ndarray_view()
        np.testing.assert_allclose(diff, 0.0, rtol = 0.0, atol = 0.0, **kwargs)
    else:
        assert_rel_equal(a1, a2, decimals, **kwargs)

def ytdataset_compare(fn1, fn2, compare_func=None, decimals = None, **kwargs):
    """
    Compare all datasets between two yt datasets.
    """

    if compare_func is None:
        compare_func = assert_array_equal

    ds1 = yt.load(fn1)
    ds2 = yt.load(fn2)

    assert ds1.field_list == ds2.field_list, \
      "Files have different datasets!"

    _kwargs = copy.deepcopy(kwargs)

    for field in ds1.field_list:
        if callable(decimals):
            _kwargs['decimals'] = decimals(field)
        elif decimals is not None:
            _kwargs['decimals'] = decimals

        compare_func(
            ds1.data[field], ds2.data[field],
            err_msg=f"Comparison of {field} field failed.", **_kwargs)
