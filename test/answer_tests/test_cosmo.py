import os
import yt

from answer_testing import \
    EnzoETest, \
    ytdataset_test, \
    assert_array_rel_equal, \
    uses_grackle, \
    cached_opts

# Set test tolerance based on compile precision
if cached_opts().uses_double_prec:
    decimals = 12
else:
    decimals = 6

# The following test is somewhat broad, but its much better than nothing. Going
# forward, we should introduce some additional tests.
    
@uses_grackle
class TestCosmoDDMultiSpecies(EnzoETest):
    parameter_file = "Cosmology/test_cosmo-dd-multispecies-short.in"
    # Note on 3/7/23: once we finish refactoring the way that grackle
    # parameters are parsed and PUPped within Enzo-E (to make use of grackle's
    # new dynamic-api), we should declare a new gold standard git-tag and use
    # parallelism in this test. Currently, if we try to run an older version of
    # Enzo-E in parallel that is linked against a newer version of Grackle the
    # PUP routines won't properly communicate Grackle parameters between the
    # different PEs
    max_runtime = 300
    ncpus = 1

    @ytdataset_test(assert_array_rel_equal, decimals=decimals)
    def test_cosmo_dd_multispecies(self):
        ds = yt.load("Dir_COSMO_MULTI_0115/Dir_COSMO_MULTI_0115.block_list")
        ad = ds.all_data()

        fluid_wfield = ("gas", "mass")
        particle_wfield = ("all", "particle_mass")
        data = {'current_time' : ds.current_time}
        for field in ds.field_list:
            if field[0] == 'enzoe':
                key, wfield = f'fluid-{field[1]}', fluid_wfield
            elif field[0] == 'all':
                key, wfield = f'particle-{field[1]}', particle_wfield
            else:
                continue
            data[key] = ad.quantities.weighted_standard_deviation(field, wfield)
        return data
