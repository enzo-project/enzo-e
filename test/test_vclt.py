import yt

from answer_testing import \
    EnzoETest, \
    ytdataset_test, \
    assert_array_rel_equal

class TestHLLCCloud(EnzoETest):
    parameter_file = "vlct/dual_energy_cloud/hllc_cloud.in"
    max_runtime = 30
    ncpus = 1

    @ytdataset_test(assert_array_rel_equal, decimals=13)
    def test_hllc_cloud(self):
        ds = yt.load("hllc_cloud_0.0625/hllc_cloud_0.0625.block_list")
        ad = ds.all_data()

        wfield = ("gas", "mass")
        data = {field[1]: ad.quantities.weighted_standard_deviation(field, wfield)
                for field in ds.field_list}

        return data
