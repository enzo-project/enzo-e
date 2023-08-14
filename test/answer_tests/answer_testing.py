import copy
from dataclasses import dataclass
import numpy as np
import os
import pytest
import shutil
import tempfile
from typing import Optional
import yt

from numpy.testing import assert_array_equal
from unittest import TestCase
from yt.funcs import ensure_dir
from yt.testing import assert_rel_equal

from test_utils.enzoe_driver import EnzoEDriver

_base_file = os.path.basename(__file__)

@dataclass(frozen = True)
class TestOptions:
    enzoe_driver: EnzoEDriver
    uses_double_prec: bool
    generate_results: bool
    test_results_dir: str
    grackle_input_data_dir : Optional[str]

_CACHED_OPTS = None

def set_cached_opts(**kwargs):
    global _CACHED_OPTS
    if _CACHED_OPTS is not None:
        raise RuntimeError("Can't call set_cached_opts more than once")

    _CACHED_OPTS = TestOptions(**kwargs)
    yt.mylog.info(
        f"{_base_file}: generate_results = {_CACHED_OPTS.generate_results}")

    yt.mylog.info(
        f"{_base_file}: test_results_dir = {_CACHED_OPTS.test_results_dir}")
    if _CACHED_OPTS.generate_results:
        ensure_dir(_CACHED_OPTS.test_results_dir)
    elif not os.path.exists(_CACHED_OPTS.test_results_dir):
        raise RuntimeError(
            f"Test results dir not found: {_CACHED_OPTS.test_results_dir}.")

    if ((_CACHED_OPTS.grackle_input_data_dir is not None) and
        (not os.path.exists(_CACHED_OPTS.grackle_input_data_dir))):
        raise RuntimeError(
            "grackle input data dir not found: "
            f"{_CACHED_OPTS.grackle_input_data_dir}")

    yt.mylog.info(
        f"{_base_file}: use_double = {_CACHED_OPTS.uses_double_prec}")

def cached_opts():
    if _CACHED_OPTS is None:
        raise RuntimeError("set_cached_opts was never called")
    return _CACHED_OPTS

src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))

input_dir = "input"

_grackle_tagged_tests = set()

def uses_grackle(cls):
    """
    Decorator that annotates that a test class uses grackle

    In detail, this annotates and sets up the appropriate skipif marker
    and updates a global registry of tests using grackle.
    """
    _grackle_tagged_tests.add(cls.__name__)

    has_grackle = cached_opts().enzoe_driver.query_has_grackle()
    has_grackle_inputs = cached_opts().grackle_input_data_dir is not None

    skip_reason = "Enzo-E is not built with Grackle"
    if has_grackle and (not has_grackle_inputs):
        skip_reason = "the grackle input data dir was not specified"

    wrapper_factory = pytest.mark.skipif(
        (not has_grackle) or (not has_grackle_inputs), reason = skip_reason
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
            with os.scandir(cached_opts().grackle_input_data_dir) as it:
                for entry in it:
                    if not entry.is_file():
                        continue
                    os.symlink(entry.path,os.path.join(self.tmpdir, entry.name))

    def setUp(self):
        self.curdir = os.getcwd()
        self.tmpdir = tempfile.mkdtemp()

        self.setup_symlinks()
        os.chdir(self.tmpdir)
        cached_opts().enzoe_driver.run(
            parameter_fname = os.path.join(input_dir, self.parameter_file),
            max_runtime = self.max_runtime, ncpus = self.ncpus,
            sim_name = f"Simulation {self.__class__.__name__}")

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
            result_filename = os.path.join(cached_opts().test_results_dir,
                                           filename)

            # check that answers exist
            if not cached_opts().generate_results:
                assert os.path.exists(result_filename), \
                  "Result file, %s, not found!" % result_filename

            data = func(*args)
            fn = yt.save_as_dataset({}, filename=filename, data=data)

            # if generating, move files to results dir
            if cached_opts().generate_results:
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
