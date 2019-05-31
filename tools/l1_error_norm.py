import argparse
import os.path
import numpy as np
import h5py
from load_enzoe import load_enzoe,load_full_grid

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

def _sanitize_units(arr):
    if isinstance(arr,yt.YTArray):
        # presumably they are already in code units
        return arr.value
    return arr

def L1_error_vector(ds,ds2,Nx=None,Ny=None,Nz=None,
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
    
    for key in ["density","velocity_x","velocity_y","velocity_z","pressure",
                "bfield_x","bfield_y","bfield_z"]:
        
        if Nx is None:
            normalize = ad[key].size
        else:
            normalize = Nx*Ny*Nz

        temp = np.abs(_sanitize_units(data[key]) - 
                      _sanitize_units(data2[key]))

        residuals.append((np.sum(np.abs(_sanitize_units(data[key]) - 
                                        _sanitize_units(data2[key])),
                                 axis = sum_axis)/float(normalize)))

    return np.array(residuals)

def norm_L1_error(ds,ds2,Nx=None,Ny=None,Nz=None,
                 sum_axis = None):
    # sum_axis controls which spatial axis we sum along
    err_vector = L1_error_vector(ds, ds2, Nx, Ny, Nz, sum_axis)
    return np.sqrt(np.sum(np.square(err_vector),axis=0))

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

def compare_to_1D_reference(ds, tab_fname, problem_ax):

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

    
    if 'time' in meta_data:
        # do a time comparison
        assert ds.current_time == float(meta_data["time"])
        
    if "left_edge" in meta_data and "right_edge" in meta_data:
        left_edge = float(meta_data["left_edge"])
        right_edge = float(meta_data["right_edge"])
    else:
        raise ValueError("Reference Table must include a line "
                         "with the comment left edge and right edge")

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

    # get the fixed resolution grid    
    ad = build_3D_grid(ds,sub_grid_edges)

    
    # Print total L1 Norm
    print(norm_L1_error(ad,ref))
    
    # compute L1_norms of axes in parallel
    L1_norms = norm_L1_error(ad,ref,sum_axis = ax_ind,
                             **kwargs)
    
    
    print("Std Dev of parallel axes: {}".format(np.std(L1_norms)))

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
    
    print(norm_L1_error(ds1,ds2,intersect,dim_length, dim_length, dim_length))
