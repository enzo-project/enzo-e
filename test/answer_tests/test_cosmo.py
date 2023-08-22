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


def _load_zeldovich_pancake_test_data(dirname, axis = 'x'):
    # helper function to assist with analyzing Zel'dovich pancake test problem
    #
    # This is currently very crude:
    # - we effectively only track a 1D problem that we run multiple times in
    #   parallel (Each one is run along the x-axis, but at a different y,z
    #   position).
    # - at the moment, we are just taking one of these 1D problems. But, maybe
    #   we should look at the mean and standard-deviations of this problem?

    if axis == 'x':
        idx = (slice(0,None), 0, 0)
    elif axis == 'y':
        idx = (0, slice(0,None), 0)
    elif axis == 'z':
        idx = (0, 0, slice(0,None))
    else:
        raise RuntimeError("problem with call to helper function")

    ds = yt.load(f'./{dirname}/{dirname}.block_list')
    grid = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)

    # SANITY CHECK ABOUT FINAL REDSHIFT:
    if ds.current_redshift > 0.001:
        raise RuntimeError(
            "Something is wrong! Redshift should be 0, not "
            f"{ds.current_redshift}"
        )

    return {
        'current_time' : ds.current_time,
        ('enzoe', 'density') : grid['enzoe','density'][idx],
        ('enzoe', 'velocity_x') : grid['enzoe','velocity_x'][idx],
        ('enzoe', 'internal_energy') : grid['enzoe','internal_energy'][idx],
        ('enzoe', 'total_energy') : grid['enzoe','total_energy'][idx],
    }
        

class TestZeldovichPancakeVL(EnzoETest):
    # it would be nice to make this test a little more rigorous
    # - see _load_zeldovich_pancake_test_data for some ideas
    # - we currently only test a problem along the x-axis. What if we run it
    #   along the y-axis or z-axis?

    parameter_file = "Cosmology/ZeldovichPancake/vl-HD-x-aligned.in"

    max_runtime = 20
    ncpus = 1

    @ytdataset_test(assert_array_rel_equal, decimals=decimals)
    def test_x_zeldovich_pancake_vl(self):
        return _load_zeldovich_pancake_test_data(
            'method-vl-HD_x-aligned-pancake_7.860e+01', axis = 'x'
        )

class TestZeldovichPancakePPM(EnzoETest):
    # it would be nice to make this test a little more rigorous
    # - see _load_zeldovich_pancake_test_data for some ideas
    # - we currently only test a problem along the x-axis. What if we run it
    #   along the y-axis or z-axis?

    parameter_file = "Cosmology/ZeldovichPancake/ppm-HD-x-aligned.in"

    max_runtime = 20
    ncpus = 1

    @ytdataset_test(assert_array_rel_equal, decimals=decimals)
    def test_x_zeldovich_pancake_ppm(self):
        return _load_zeldovich_pancake_test_data(
            'method-ppm-HD_x-aligned-pancake_7.860e+01', axis = 'x'
        )
