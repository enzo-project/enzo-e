import numpy as np
import os
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
yt.mylog.info(f"{_base_file}: {generate_results=}")

_results_dir = os.environ.get("TEST_RESULTS_DIR", "~/enzoe_test_results")
test_results_dir = os.path.abspath(os.path.expanduser(_results_dir))
yt.mylog.info(f"{_base_file}: {test_results_dir=}")
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
yt.mylog.info(f"{_base_file}: {charmrun_path=}")
if not os.path.exists(charmrun_path):
    raise RuntimeError(
        f"No charmrun executable found in {_charm_path}.")

src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
enzo_path = os.path.join(src_path, "build/bin/enzo-e")
yt.mylog.info(f"{_base_file}: {enzo_path=}")
if not os.path.exists(enzo_path):
    raise RuntimeError(
        f"No enzo-e executable found in {enzo_path}.")

input_dir = "input"

class EnzoETest(TestCase):
    parameter_file = None
    max_runtime = np.inf
    ncpus = None

    def setup_symlinks(self):
        ipath = os.path.join(src_path, input_dir)
        spath = os.path.join(self.tmpdir, input_dir)
        os.symlink(ipath, spath)

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
                os.rename(filename, result_filename)
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
    """
    assert_rel_equal(a1, a2, decimals, **kwargs)

def ytdataset_compare(fn1, fn2, compare_func=None, **kwargs):
    """
    Compare all datasets between two yt datasets.
    """

    if compare_func is None:
        compare_func = assert_array_equal

    ds1 = yt.load(fn1)
    ds2 = yt.load(fn2)

    assert ds1.field_list == ds2.field_list, \
      "Files have different datasets!"

    for field in ds1.field_list:
        compare_func(
            ds1.data[field], ds2.data[field],
            err_msg=f"Comparison of {field} field failed.", **kwargs)
