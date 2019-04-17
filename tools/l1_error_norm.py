import argparse
import os.path
import numpy as np
from load_enzoe import load_enzoe

parser = argparse.ArgumentParser(description='Computes norm of the L1 error ')
parser.add_argument('data_1', action = 'store', 
                    help = ("The directory containing the data that serves as "
                            "reference values for the calculation"))
parser.add_argument('data_2', action = 'store', 
                    help = ("The directory containing the data which we are "
                            "compute the L1 norm."
                            "reference values for the calculation"))
parser.add_argument('-v', action = 'store_true', default = False,
                    help = ("Indicates that messages should be verbose and "
                            "include utilized fields."))
parser.add_argument('-n', action = 'store', default = None, required = False,
                    type = int,
                    help = ("By default the underlying arrray shape is used "
                            "for normalization of the L1 norm. If a positive "
                            "integer is specified for this argument, then the "
                            "normalization is calculated using the value as "
                            "the length of each dimension."))

# should allow for specification of fields

def get_block_list(dir_name):
    if os.path.basename(dir_name) == '':
        # make sure that there are no slashes at the end of the dir_name
        stop_slice = len(dir_name)
        for i in range(stop_slice-1,-1,-1):
            if dir_name[i] == '/':
                stop_slice-=1
            else:
                break
        dir_name = dir_name[:stop_slice]
            
            
    fname = os.path.join(dir_name,
                         '.'.join((os.path.basename(dir_name), 'block_list')))
    assert os.path.isfile(fname)
    return fname


def norm_L1_error(ds,ds2,field_list,dim_length = None):
    ad = ds.all_data()
    ad2 = ds2.all_data()
    residuals = []
    for key in field_list:
        if dim_length is None:
            norm = ad[key].size
        else:
            norm = dim_length ** ds.dimensionality
        residuals.append((np.sum(np.abs(ad[key] - ad2[key]).value)/
                          float(norm)))
    return np.sqrt(np.sum(np.square(residuals)))



if __name__ == '__main__':
    
    args = parser.parse_args()

    verbose = args.v

    if not verbose:
        from yt.funcs import mylog
        mylog.setLevel(30)

    file_1 = get_block_list(args.data_1)
    file_2 = get_block_list(args.data_2)

    dim_length = args.n
    if dim_length is not None:
        assert dim_length > 0, "-n must point to a positive int"

    ds1 = load_enzoe(file_1)
    ds2 = load_enzoe(file_2)

    intersect = [elem for elem in ds1.field_list if elem in ds2.field_list]

    if verbose:
        print("field_list: {:}".format(str(intersect)))
    print(norm_L1_error(ds1,ds2,intersect,dim_length = dim_length))
