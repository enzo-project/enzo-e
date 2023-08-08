import os
import numpy as np
import yt

from answer_testing import \
    EnzoETest, \
    ytdataset_test, \
    assert_array_rel_equal, \
    uses_grackle, \
    cached_opts

# Set test tolerance based on compile precision
use_double = cached_opts().uses_double_prec
if use_double:
    decimals = 12
else:
    decimals = 6


def _get_decimals(field_name):
    # returns infinity (require an exact match), when we are looking at the
    # location of a minimum/maximum.
    assert isinstance(field_name, tuple)
    if field_name[1].split(':')[-1] not in ["min", "max", "mean", "std_dev"]:
        return np.inf
    elif use_double:
        return 14
    else:
        return 6

@uses_grackle
class TestGrackleGeneral(EnzoETest):
    parameter_file = "Grackle/method_grackle_general.in"
    # Note on 3/7/23: once we finish refactoring the way that grackle
    # parameters are parsed and PUPped within Enzo-E (to make use of grackle's
    # new dynamic-api), we should declare a new gold standard git-tag and use
    # parallelism in this test. Currently, if we try to run an older version of
    # Enzo-E in parallel that is linked against a newer version of Grackle the
    # PUP routines won't properly communicate Grackle parameters between the
    # different PEs
    max_runtime = 240
    ncpus = 1

    @ytdataset_test(assert_array_rel_equal, decimals = _get_decimals)
    def test_grackle_general(self):
        """
        Compares a variety of summary statistics.

        This was adapted from an earlier test.
        """
        ds = yt.load("GeneralGrackle-500.00/GeneralGrackle-500.00.block_list")
        ad = ds.all_data()

        quan_entry_sets = [
            ('min_location', ('min', 'min_xloc', 'min_yloc', 'min_zloc'), {}),
            ('max_location', ('max', 'max_xloc', 'max_yloc', 'max_zloc'), {}),
            ('weighted_standard_deviation',
             ('std_dev', 'mean'), {'weight' : ("gas", "cell_volume")}),
        ]

        data = {}
        for field in ds.field_list:
            for derived_quantity, rslt_names, kwargs in quan_entry_sets:
                quan_func = getattr(ad.quantities,derived_quantity)
                rslts = quan_func(field, **kwargs)
                num_rslts = len(rslts)
                assert num_rslts == len(rslt_names)
                for i in range(num_rslts):
                    data[f'{field[1]}:{rslt_names[i]}'] = rslts[i]
        return data

@uses_grackle
class TestGrackleCoolingDt(EnzoETest):
    parameter_file = "Grackle/method_grackle_cooling_dt.in"
    # Note on 3/7/23: once we finish refactoring the way that grackle
    # parameters are parsed and PUPped within Enzo-E (to make use of grackle's
    # new dynamic-api), we should declare a new gold standard git-tag and use
    # parallelism in this test. Currently, if we try to run an older version of
    # Enzo-E in parallel that is linked against a newer version of Grackle the
    # PUP routines won't properly communicate Grackle parameters between the
    # different PEs
    max_runtime = 120
    ncpus = 1

    @ytdataset_test(assert_array_rel_equal, decimals=decimals)
    def test_grackle_cooling_dt(self):
        """
        Compares the current time at cycle 20.

        This was adapted from an earlier test that just compared the completion
        time of the test. The functionallity to just check completion times did
        not exist while adding this test, so we need to read in a full output.

        It may be worth introducing this functionallity in the future (it could
        be acheived by simply parsing the logs). Alternatively, it might be
        worth reworking this test.
        """
        ds = yt.load(
            "GrackleCoolingDt-cycle20/GrackleCoolingDt-cycle20.block_list")
        return {'current_time' : ds.current_time}
