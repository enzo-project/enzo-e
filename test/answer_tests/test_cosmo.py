import os
import yt

from answer_testing import \
    EnzoETest, \
    ytdataset_test, \
    assert_array_rel_equal, \
    uses_grackle

_base_file = os.path.basename(__file__)

# Set test tolerance based on compile precision
use_double = os.environ.get("USE_DOUBLE", "false").lower() == "true"
yt.mylog.info(f"{_base_file}: use_double = {use_double}")
if use_double:
    decimals = 12
else:
    decimals = 6

# The following test is somewhat broad, but its much better than nothing. Going
# forward, we should introduce some additional tests.
    
@uses_grackle
class TestCosmoDDMultiSpecies(EnzoETest):
    parameter_file = "Cosmology/test_cosmo-dd-multispecies-short.in"
    max_runtime = 75
    ncpus = 4

    @ytdataset_test(assert_array_rel_equal, decimals=decimals)
    def test_cosmo_dd_multispecies(self):
        ds = yt.load("Dir_COSMO_MULTI_0115/Dir_COSMO_MULTI_0115.block_list")
        ad = ds.all_data()

        # TODO: refactor this code so that we compare standard deviations
        # - at the time of writing this, there seems to be a bug with invoking
        #   ad.quantities_standard_deviation on a simulation with AMR...
        #wfield = ("gas", "mass")
        #data = {
        #    field[1]: ad.quantities.weighted_standard_deviation(field, wfield)
        #    for field in ds.field_list}
        return {'current_time' : ds.current_time}
