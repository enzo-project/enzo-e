import argparse
import os.path
import warnings
import numpy as np

warnings.simplefilter("ignore", FutureWarning)
import yt
import h5py
from load_enzoe import load_enzoe,load_full_grid
warnings.resetwarnings()


_description = '''\
Calculates the norm of the L1 error vector of enzo-e results.

Presently, this program only supports 3D unigrid enzo-e data and it runs in 
2 modes: "sim" and "table" (reflecting different types of reference data). See 
the end for more details.
'''

_epilog = '''\
sim mode
=========
  In "sim" mode the reference data comes from another enzo-e dataset with the 
  same layout as the target data (This might be useful for perodic wave 
  evolution). This calculation is inherently 3D'

  Example Usage:
    l1_error_norm.py path/to/ref/data path/to/target/data

table mode
==========
  In "table" mode the reference data comes from a table stored in a csv file. 
  Presently this mode assumes that problem is inherently 1D along a given axis 
  (The table should holds data that is expected for all possible 1D slices 
  along the given axis). This axis must be specified as x, y, or z immediately
  after the mode is indicated. 

  By default, the yielded L1 error norm uses the data from all 3D data in the 
  active region of the problem. Specifying the --std option will yield the 
  standard deviation of the L1 error norms computed for every 1D slice along 
  the active axis

  Example Usage:
    Total L1 Error Norm for a Riemann problem evolved along the y-axis:\n'
    $ python l1_error_norm.py table y path/to/table path/to/target/data

    Standard Deviation of L1 Error Norm for all 1D y-slices
    $ python l1_error_norm.py table y --std path/to/table path/to/target/data

  Table Format: 
    The supplied table should be comma separated and header comments must be led
    by a "#". 2 separate lines in the header file must be dedicated to 
    specifying
      #left_edge = ...
      #right_edge = ...
    to specify the the left and right grid boundaries along the active axis 
    where the table values should be compared (The domain of the grid may 
    extend further along the active region). (These comments can be specified 
    in any order and can be mixed with other comments - additionally whitespace 
    is not a concern). 

    The first line following the comments should include the names of the 
    columns (they should match the fluid fluid names in the simulation to 
    compare against). There must be a single column titled position, pos, x, y, 
    or z indicating the (cell-centered) position of each value and this column 
    must increase in value monotonically.
'''


_formatter_class = argparse.RawDescriptionHelpFormatter

parser = argparse.ArgumentParser(description=_description,
                                 formatter_class = _formatter_class,
                                 epilog = _epilog)

parser.add_argument('ref_type',
                    choices = ['sim','tablex', 'tabley', 'tablez'],
                    help = ('Indicates the type of reference data being. '
                            '"sim" refers to 3D simulation data while "table" '
                            'refers to 1D table data. In the latter case, the '
                            'dimension of the target data along which the '
                            'table can be compared must also be specified (as '
                            'x, y, or z).'))


parser.add_argument('ref_path', action = 'store', 
                    help = ("In \"sim\" mode this is the directory containing "
                            "the simulation data that serves as the reference "
                            "values. In \"table\" mode this file name of the "
                            "table containing the reference values"))
parser.add_argument('target_path', action = 'store', 
                    help = ("The directory containing the simulation data for "
                            "which the L1 Error Norm is computed."))
parser.add_argument('-v', action = 'store_true', default = False,
                    help = ("Indicates that messages should be verbose and "
                            "include utilized fields."))
parser.add_argument('-n', action = 'store', default = None, required = False,
                    type = int,
                    help = ("This only applies to \"sim\" mode. By default the "
                            "underlying grid shape is used for normalization "
                            "of the L1 norm. If a positive integer is "
                            "specified for this option, then the normalization "
                            "uses this value as the length of each dimension."))
parser.add_argument('--std', action = 'store_true', default = False,
                    help = ("This only applies to \"table\" mode. By default, "
                            "L1 error norm is computed for the 3D grid. If "
                            "this option is supplied, then the standard "
                            "deviation of the L1 error norms for all 1D slices "
                            "is returned."))

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

def _sanitize_units(arr):
    if isinstance(arr,yt.YTArray):
        # presumably they are already in code units
        return arr.value
    return arr

def L1_error_vector(ds, ds2, comparison_fields, Nx=None, Ny=None, Nz=None,
                    sum_axis = None):
    # sum_axis indicate the axis along which to sum along

    assert ((Nx is None and Ny is None and Nz is None) or 
            (Nx is not None and Ny is not None and Nz is not None))

    if isinstance(ds,yt.data_objects.static_output.Dataset):
        data = ds.all_data()
    else:
        data = ds

    if isinstance(ds2,yt.data_objects.static_output.Dataset):
        data = ds2.all_data()
    else:
        data2 = ds2

    residuals = []

    for key in comparison_fields:
        
        if Nx is None:
            normalize = data[key].size
        else:
            normalize = Nx*Ny*Nz

        temp = np.abs(_sanitize_units(data[key]) - 
                      _sanitize_units(data2[key]))

        residuals.append((np.sum(np.abs(_sanitize_units(data[key]) - 
                                        _sanitize_units(data2[key])),
                                 axis = sum_axis)/float(normalize)))

    return np.array(residuals)

def norm_L1_error(ds, ds2, comparison_fields, Nx=None, Ny=None, Nz=None,
                 sum_axis = None):
    # sum_axis controls which spatial axis we sum along
    err_vector = L1_error_vector(ds, ds2, comparison_fields, Nx, Ny, Nz,
                                 sum_axis)
    return np.sqrt(np.sum(np.square(err_vector),axis=0))

def find_common_fields(ds1,ds2, verbose):
    if isinstance(ds1,yt.data_objects.static_output.Dataset):
        field_list1 = [elem[1] for elem in ds1.field_list]
    else:
        field_list1 = ds1.keys()

    if isinstance(ds2,yt.data_objects.static_output.Dataset):
        field_list2 = [elem[1] for elem in ds2.field_list]
    else:
        field_list2 = ds2.keys()

    intersect = [elem for elem in field_list1 if elem in field_list2]

    if verbose:
        print("field_list: {:}".format(str(intersect)))
    return intersect


# Functions used to compute the L1Error Norm between a simulation and 1D
# reference table

# Absolute tolerance for the error in the specified cell boundary location
_ATOL = 1.e-13

def _find_boundary_index(val, cell_width, domain_left_edge,ax_num):
    """
    Finds the index of the cell who's left boundary 
    corresponds to val.
    """
    if isinstance(val,yt.YTQuantity):
        start = val.to('code_length').value
    else:
        start = float(val)
        
    index = int((start-domain_left_edge.value) /cell_width.value + 0.5)
    if not np.isclose((domain_left_edge+index*cell_width).value,
                       val, atol=_ATOL, rtol=0):
        raise ValueError("{} does not lie at a cell boundary".format(val)
                         + " along axis {}".format(ax_num))
    return index
    

def build_3D_grid(ds, dim_edges):
    """
    Builds a CoveringGrid
    
    Parameters
    ----------
    ds: ~yt.data_objects.static_output.Dataset
        The dataset to make the grid from
    dim_edges : sequence
        A sequence of 3 entries where each entry indicates what 
        portion of an axis to include. An entry should either be None,
        indicating that we should include the full dimension or a list 
        of the start and stop positions (which must line up with cell 
        boundaries)
    """
    # reversed_order applies to custom load_enzoe function
    if len(dim_edges) != 3:
        raise ValueError("dim_edges must have 3 elements")
    grid_dim = ds.domain_dimensions
    
    cell_widths = ds.domain_width/ds.domain_dimensions
    
    iter_tuple = zip(dim_edges, ds.domain_left_edge,
                     ds.domain_dimensions, cell_widths)
    
    grid_left_edge = []
    dims = []
    for i,temp in enumerate(iter_tuple):
        elem, domain_left_edge, ax_length, width = temp
        
        if elem is None:
            grid_left_edge.append(domain_left_edge)
            dims.append(ax_length)
            continue
        assert(len(elem)==2)
        start_ind = _find_boundary_index(elem[0], width, 
                                         domain_left_edge, i)
        grid_left_edge.append(elem[0])
        stop_ind = _find_boundary_index(elem[1], width, 
                                         domain_left_edge, i)
        dims.append(stop_ind-start_ind)
    return ds.covering_grid(0,left_edge = grid_left_edge, dims = dims)
        

def load_reference_table(fname):
    
    skip = 0
    meta_data = {}
    with open(fname,"r") as f:
        for line in f.readlines():
            if (line[0]) == "#":
                skip+=1
                if "=" in line:
                    words = line[1:].split("=")
                    if len(words) == 2:
                        meta_data[words[0].strip()] = words[1].strip()
            else:
                break
    rec = np.recfromtxt(fname, skip_header=skip,comments="#",dtype=None,
                        delimiter=',',names=True,encoding='utf-8')
    
    out_data = {}
    for field in rec.dtype.fields:
        out_data[field] = rec[field]
    return meta_data, out_data

def compare_to_1D_reference(ds, tab_fname, problem_ax, verbose, calc_std):

    meta_data, ref_data = load_reference_table(tab_fname)

    pos_keys = ["position","pos", "x", "y", "z"]
    intersect = [val for val in ref_data.keys() if val in pos_keys]
    if len(intersect) == 0:
        raise ValueError("table does not have a column labelled "
                         "position, pos, x, y, or z")
    elif len(intersect) >1:
        raise ValueError("table must have only one column labelled "
                         "position, pos, x, y, or z")
    else:
        pos = ref_data.pop(intersect[0])
        ref_data["pos"] = pos

    #left_edge = (pos[0] - np.diff(pos)[0]*0.5)
    #right_edge = (pos[-1] + np.diff(pos)[0]*0.5)

    if 'time' in meta_data:
        # do a time comparison
        assert ds.current_time == float(meta_data["time"])
        
    if "left_edge" in meta_data and "right_edge" in meta_data:
        left_edge = float(meta_data["left_edge"])
        right_edge = float(meta_data["right_edge"])
    else:
        raise ValueError("Reference Table must include a line "
                         "with the comment left edge and right edge")

    #kwargs is only used if computing standard deviation
    kwargs = {"Nx":1, "Ny":1, "Nz":1}
    if problem_ax == "z":
        ax_ind = 2
        kwargs["Nz"] = len(pos)
    elif problem_ax == "y":
        ax_ind = 1
        kwargs["Ny"] = len(pos)
    else:
        ax_ind = 0
        kwargs["Nx"] = len(pos)

    # adjust the shape of the reference arrays
    # and determine limits along problem_ax that we are interested in
    temp = []
    sub_grid_edges = []
    for i in range(3):
        if i != ax_ind:
            temp.append(np.newaxis)
            sub_grid_edges.append(None)
        else:
            temp.append(slice(None,None))
            sub_grid_edges.append([left_edge,right_edge])
    idx = tuple(temp)
    ref = {}
    for field,arr in ref_data.items():
        ref[field] = arr[idx]  

    # get the fixed resolution grid of the target dataset
    ad = build_3D_grid(ds,sub_grid_edges)

    comparison_fields = find_common_fields(ds,ref, verbose)

    if calc_std:
        # compute L1_norms of axes in parallel
        L1_norms = norm_L1_error(ad, ref, comparison_fields, sum_axis = ax_ind,
                                 **kwargs)
        print(np.std(L1_norms))
    else:
        # just compute the L1 Error Norm normally
        print(norm_L1_error(ad,ref, comparison_fields))

# Functions used to compute L1Error Norm between simulations

def compare_to_sim_reference(ds, ref_path, verbose, dim_length):
    ds_ref = load_enzoe(get_block_list(ref_path))
    comparison_fields = find_common_fields(ds_ref, ds, verbose)
    print(norm_L1_error(ds_ref, ds, comparison_fields, Nx = dim_length,
                        Ny = dim_length, Nz=dim_length))

if __name__ == '__main__':
    args = parser.parse_args()

    verbose = args.v
    if not verbose:
        from yt.funcs import mylog
        mylog.setLevel(30)

    # load the target data to be analyzed
    target_file = get_block_list(args.target_path)
    ds = load_enzoe(target_file)

    if args.ref_type == 'sim':
        sim_comp = True
        if args.std:
            raise ArgumentError('The --std option cannot be specified for '
                                '"sim" reference data.')
        dim_length = args.n
        if (dim_length is not None) and (dim_length <= 0):
            raise ArgumentError('The -n option must be followed by a positive '
                                'integer.')
        compare_to_sim_reference(ds, args.ref_path, verbose, dim_length)

    else:
        sim_comp = False
        if args.n is not None:
            raise ArgumentError('The -n option cannot be specified for '
                                '"table" reference data')
        
        problem_ax = args.ref_type[-1]
        compare_to_1D_reference(ds, args.ref_path, problem_ax, verbose,
                                args.std)
