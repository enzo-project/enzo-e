import argparse
import os.path
import warnings
import numpy as np

warnings.simplefilter("ignore", FutureWarning)
import yt
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

  By default, the yielded L1 error norm for the all of the 3D data in the 
  problem region (this is equivalent to taking the average of the L1-Norms of
  for every 1D slice along the active axis. To take the median or standard 
  deviations of these L1-Norms pass the --median OR --std options

  Example Usage:
    Total L1 Error Norm for a Riemann problem evolved along the y-axis:\n'
    $ python l1_error_norm.py tabley path/to/table path/to/target/data

    Standard Deviation of L1 Error Norm for all 1D y-slices
    $ python l1_error_norm.py tabley --std path/to/table path/to/target/data

  Table Format:
    The supplied table should be comma separated and header comments must be led
    by a "#". The first line following the comments should include the names of 
    the columns (they should match the fluid fluid names in the simulation to 
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
                    help = ('Indicates the type of reference data being used. '
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

class StoreCommaListAction(argparse.Action):
    """Stores comma a list from comma separated values (no spaces). 

    It's okay if there is a trailing comma. A list of the strings are stored
    at the destination.
    """
    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        if nargs is not None:
            raise ValueError("nargs not allowed")
        super(StoreCommaListAction, self).__init__(option_strings, dest,
                                                   **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        if values[-1] == ',':
            parsed_vals = values[:-1].split(',')
        else:
            parsed_vals = values.split(',')
        setattr(namespace, self.dest, parsed_vals)

parser.add_argument('-f','--fields', action = StoreCommaListAction,
                    default = None, required = False,
                    help = ("Specify the list of fields to use to compute the "
                            "error norm. The field names should be separated "
                            "by commas and have no spaces."))
parser.add_argument('-n', action = 'store', default = None, required = False,
                    type = int,
                    help = ("This only applies to \"sim\" mode. By default the "
                            "underlying grid shape is used for normalization "
                            "of the L1 norm. If a positive integer is "
                            "specified for this option, then the normalization "
                            "uses this value as the length of each dimension."))
parser.add_argument('--permute', nargs = '?', action = 'store', default = None,
                    const = '-1', type = int,
                    help = ("This only applies to \"table\" mode. By default "
                            "the column names from the reference table are "
                            "assumed to match the compared simulation field "
                            "names. A value of 1 or 2 will cyclically permute "
                            "the axes of vector quantities (any collection of "
                            "column names with a common prefix and suffixes of "
                            "\"_x\", \"_y\" and \"_z\"). If the option is "
                            "passed but no argument is used, a guess will be "
                            "made based on the specified active axis"))
parser.add_argument('--std', action = 'store_true', default = False,
                    help = ("This only applies to \"table\" mode. By default, "
                            "L1 error norm is computed for the 3D grid. If "
                            "this option is supplied, then the standard "
                            "deviation of the L1 error norms for all 1D slices "
                            "is returned."))
parser.add_argument('--median', action = 'store_true', default = False,
                    help = ("This only applies to \"table\" mode. By default, "
                            "L1 error norm is computed for the 3D grid. If "
                            "this option is supplied, then the median of the "
                            "L1 error norms for all 1D slices is returned."))

parser.add_argument('--reverse', action = 'store_true', default = False,
                    help = ("This only applies to \"table\" mode. When "
                            "specified, the position values in the table are "
                            "reversed and all vector components are "
                            "multiplied by -1."))

parser.add_argument('--offset_soln', action = 'store', default = 0.0,
                    type= float,
                    help = ("This only applies to \"table\" mode. When "
                            "specified, this value is added to the positions "
                            "specified in the table along the comparison "
                            "axis. This must be evenly divisable by the cell "
                            "width."))

parser.add_argument('--bkg_velocity', action = 'store', nargs = 3,
                    default = [], type = float,
                    help = ("This only applies to \"table\" mode. This "
                            "option is used to specify constant background "
                            "velocities that are not specified in the table."
                            "This information is used to update the expected "
                            "velocities and total energy. This option "
                            "requires 3 arguments to specify x,y,z velocity "
                            "components. This is unaffected by the --permute "
                            "and --reverse options."))

# should allow for specification of fields
def get_block_list(dir_name):
    orig_dir_name = dir_name
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
    if not os.path.isfile(fname):
        if not os.path.isdir(dir_name):
            raise ValueError(
                "{:s} is not the name of a directory".format(orig_dir_name))
        raise ValueError("a file called {:s} can't be found".format(fname))
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
        data2 = ds2.all_data()
    else:
        data2 = ds2

    residuals = []

    for key in comparison_fields:
        
        if Nx is None:
            normalize = data[key].size
        else:
            normalize = Nx*Ny*Nz

        # sorting the absolute differences before summing them would allow for
        # better symmetry tests
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

_DEFAULT_COMMON_FIELDS \
    = ["density", "total_energy", "velocity_x", "velocity_y", "velocity_z",
       "bfield_x", "bfield_y", "bfield_z"]

def find_common_fields(ds1, ds2, verbose, fields = None):
    """
    Returns a list of fields common to ds1 and ds2

    Parameters
    ----------
    ds1, ds2 : `yt.data_objects.static_output.Dataset` or `dict` of `ndarray`
        These represent the different datasets. Field names are obtained for 
        the former type by accessing the `derived_field_list` attribute. For
        the latter type, they are obtained by calling the `keys` method.
    verbose : `bool`
        If True, prints the list of common fields
    fields : `list` of strings or `None` (optional)
        If specified, the function checks that both datasets contain this list 
        of fields (this is returned by the function). If `None` (default),
        this selectively returns the subset of common fields from the 
        `_DEFAULT_COMMON_FIELDS`.
    """
    if isinstance(ds1,yt.data_objects.static_output.Dataset):
        field_list1 = [elem[1] for elem in ds1.derived_field_list]
    else:
        field_list1 = ds1.keys()

    if isinstance(ds2,yt.data_objects.static_output.Dataset):
        field_list2 = [elem[1] for elem in ds2.derived_field_list]
    else:
        field_list2 = ds2.keys()

    if fields is None:
        intersect = []
        for elem in _DEFAULT_COMMON_FIELDS:
            if elem in field_list1 and elem in field_list2:
                intersect.append(elem)
        #if len(intersect) == 0:
        #    intersect = [elem for elem in field_list1 if elem in field_list2]
    else:
        # check that both datasets have all of the user-specified fields
        for field in fields:
            if (field not in field_list1) or (field not in field_list2):
                raise ValueError("Both datasets do not have the field "
                                 +"\"{:s}\".".format(field))
        intersect = fields

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

def vector_name_iterator(data):
    """
    Produces an iterator that yields 2-tuples of vectors given a dict of fields

    vectors are identified by identifying fields with common prefixes that all
    end with '_x', '_y', or '_z'.

    The first element of the yielded tuple holds the common prefix of the 
    fields related to the vector, while the second element holds a list of 
    field names corresponding to the various components (orderred as x,y,z).
    Missing components are replaced by a None
    """

    ax_map = {'x' : 0, 'y' : 1, 'z' : 2}
    candidates = {}
    # identify candidate vectors
    for elem in data.keys():
        temp = elem.split('_')
        if len(temp) != 2 or (temp[1] not in ['x','y','z']):
            continue
        prefix, dim = temp
        if prefix not in candidates:
            candidates[prefix] = [None, None, None]
        candidates[prefix][ax_map[dim]] = elem
    return candidates.items()

def permute_vector_quantites(data, verbose, permute):

    if permute is None:
        return None
    assert permute in [0,1,2]

    vectors = []
    for vector_name, names in vector_name_iterator(data):
        if len(list(filter(lambda x: x is not None, names))) != 3:
            # skip this vector if all 3 components are not provided
            continue
        vectors.append(vector_name)
        temp = [data[names[0]], data[names[1]], data[names[2]]]
        data[names[permute]] = temp[0]
        data[names[(permute+1)%3]] = temp[1]
        data[names[(permute+2)%3]] = temp[2]

    if verbose:
        print("permuted the vectors {}".format(str(vectors)))

def reverse_1D_soln(data):
    """
    Reverses the domain of the 1D tabulated solution.
    
    Namely reverse the order of all fields and multiply the values of all 
    vector components by -1.
    """

    for key,value in data.items():
        if key != 'pos':
            data[key] = value[::-1]

    for _, vector_field_l in vector_name_iterator(data):
        for field_name in vector_field_l:
            if field_name is None:
                continue
            data[field_name] *= -1.

def compare_to_1D_reference(ds, tab_fname, problem_ax, verbose, op_func = None,
                            permute = 0, specified_fields = None,
                            reverse = False, offset_soln = 0.,
                            bkd_velocity = []):

    # if op_func is None, returns the total L1 Error Norm
    # otherwise op_func is applied on all parallel L1 Error Norms

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
        if (np.diff(ref_data[intersect[0]])<=0).any():
            msg = 'the {} column in table must monotonically increase'
            raise ValueError(msg.format(intersect[0]))
        pos = ref_data.pop(intersect[0])
        ref_data["pos"] = pos

    # if offset_soln was specified, modify the position values
    cell_width = np.diff(pos)[0]
    if offset_soln != 0.:
        if (offset_soln % np.abs(cell_width) != 0.):
            msg = ('the offset_soln value is not an integral multitple of '
                   'cell_width (which is {0!r}')
            raise ValueError(msg.format(float(cell_width)))
        ref_data["pos"] += offset_soln

    # compute left_edge and right_edge
    left_edge = (pos[0] - cell_width*0.5)
    right_edge = (pos[-1] + cell_width*0.5)

    if reverse:
        reverse_1D_soln(ref_data)

    # add in some missing vector components
    fields_to_check = ["velocity_x", "velocity_y", "velocity_z"]
    bfield_l = ['bfield_x', 'bfield_y', 'bfield_z']
    if any( [ (field_name in ref_data) for field_name in bfield_l ]):
        # only add bfield components if at least one component already existed
        fields_to_check = fields_to_check + bfield_l
    for field in fields_to_check:
        if field not in ref_data:
            ref_data[field] = np.zeros_like(pos)

    # it might be nice to compute "derived fields". Namely total_energy,
    # internal_energy, and/or pressure. We should require that gamma be
    # specified in meta_data

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

    # optionally permute some vector quantities
    permute_vector_quantites(ref, verbose, permute)

    # optionally add some fixed background velocities
    if bkd_velocity != []:
        # make sure to adjust total energy
        for i, vel_name in enumerate(['velocity_x', 'velocity_y',
                                      'velocity_z']):
            bkg_v = bkg_velocity[i]
            old_varr = ref[vel_name].copy()
            ref[vel_name] += bkg_v
            new_varr = ref[vel_name]
            if 'total_energy' in ref:
                d_ke = 0.5*(np.square(new_varr) - np.square(old_varr))
                ref['total_energy'] += d_ke

    # get the fixed resolution grid of the target dataset
    ad = build_3D_grid(ds,sub_grid_edges)

    comparison_fields = find_common_fields(ds,ref, verbose, specified_fields)

    if op_func is None:
        # just compute the L1 Error Norm normally
        print(norm_L1_error(ad, ref, comparison_fields))
    else:
        # compute L1_norms of axes in parallel
        L1_norms = norm_L1_error(ad, ref, comparison_fields, sum_axis = ax_ind,
                                 **kwargs)
        print(op_func(L1_norms))

def yt_load(fname):
    warnings.simplefilter("ignore", ResourceWarning)
    warnings.simplefilter("ignore", UserWarning)
    ds = yt.load(fname)
    warnings.resetwarnings()
    return ds

# Functions used to compute L1Error Norm between simulations
def compare_to_sim_reference(ds, ref_path, verbose, dim_length,
                             specified_fields = None):
    ds_ref = yt_load(get_block_list(ref_path))
    comparison_fields = find_common_fields(ds_ref, ds, verbose,
                                           specified_fields)
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
    ds = yt_load(target_file)

    specified_fields = args.fields

    if args.ref_type == 'sim':
        sim_comp = True

        unexpected_args = [("std",False), ("median", False), ("permute", None),
                           ("reverse",False), ("offset_soln", 0.),
                           ('bkg_velocity', [])]
        for name,default_val in unexpected_args:
            if hasattr(args,name) and getattr(args,name) != default_val:
                raise ValueError(('The --{} option cannot be specified for '
                                  '"sim" reference data.').format(name))
        dim_length = args.n
        if (dim_length is not None) and (dim_length <= 0):
            raise ValueError('The -n option must be followed by a positive '
                             'integer.')
        compare_to_sim_reference(ds, args.ref_path, verbose, dim_length,
                                 specified_fields)

    else:
        sim_comp = False
        if args.n is not None:
            raise ValueError('The -n option cannot be specified for '
                             '"table" reference data')

        op_func = None
        
        if args.std and args.median:
            raise ValueError("The --std and --median options can not both "
                             "be specified at the same time.")
        elif args.std:
            op_func = np.std
        elif args.median:
            op_func = np.median

        problem_ax = args.ref_type[-1]

        # check the argument passed to permute
        if args.permute is None:
            permute = args.permute
        elif args.permute == -1:
            # specified the option without an argument
            permute = {'x' : 0, 'y' : 1, 'z' : 2}[problem_ax]
        else:
            permute = int(args.permute)
            if permute not in [0, 1, 2]:
                raise ValueError(("argument --permute: invalid choice: {}"
                                  "(don't pass an argument or choose from "
                                  "'0', '1', '2')").format(args.permute))

        reverse = args.reverse
        offset_soln = args.offset_soln
        bkg_velocity = getattr(args,'bkg_velocity',[])
        

        compare_to_1D_reference(ds, args.ref_path, problem_ax, verbose, op_func,
                                permute, specified_fields, reverse,
                                offset_soln, bkg_velocity)
