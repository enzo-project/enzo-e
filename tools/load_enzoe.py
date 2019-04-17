
"""
This file includes the definition of load_full_grid and load_enzoe which are 
useful for reading in uniform grids Enzo-e results for regression testing.

At the time of writing, there are 2 bugs in the yt front-end that prevent us 
from just using yt for our regression tests:
    1. There are bugs with loading in enzo-e dataset when using python2 (which 
       is required for the build system)
    2. The front-end cannot currently load in the results from serial runs.

If the above 2 bugs are fixed, this file should be removed
"""
import os.path
import warnings
import numpy as np

warnings.simplefilter("ignore", FutureWarning)
import h5py
import yt
warnings.resetwarnings()

def _get_slices(attrs):
    return (slice(attrs[u'enzo_GridStartIndex'][2],
                  attrs[u'enzo_GridEndIndex'][2]+1),
            slice(attrs[u'enzo_GridStartIndex'][1],
                  attrs[u'enzo_GridEndIndex'][1]+1),
            slice(attrs[u'enzo_GridStartIndex'][0],
                  attrs[u'enzo_GridEndIndex'][0]+1))

def _prepare_array_of_blocks(block_data):
    
    # prepares a list of lists, where each entry 
    # corresponds to a tuple of the block_name and the 
    # file that contains it
    
    ndims = len(block_data[0][0][1:].split("_"))
    assert ndims == 3, "Not currently equipped to handle fewer dims"

    root_shape = [0 for elem in range(ndims)]
    
    # determine the root shape 
    
    index_l = []
    for block, fname in block_data:
        index = []
        for i,val in enumerate(block[1:].split("_")):
            if val == '':
                index_val =0
            else:
                index_val = int(val,2)
            root_shape[i] = max(root_shape[i],index_val+1)
            index.append(index_val)
        index_l.append(tuple(index))

    # now we have root_shape ordered x,y,z
    # index_l also lists the indices for an array, but list them as (x,y,z)
    # instantiate array ordered z,y,x
    block_array = np.empty(shape=root_shape[::-1],dtype = object)

    # now fill in the arrays with the appropriate tuples
    for block_tup, index in zip(block_data,index_l):
        block_array[index[::-1]] = block_tup
    
    return root_shape, block_array

def get_details_(dir_name, block_data):
    # gets basic data about the simulation
    #    - resolution (all block have the same resolution)
    #    - list of field tuples (first entry indicates name in hdf5 file,
    #      second entry indicates the final name)
    #    - left edge and right edge of grid
    #    - size of ghost zone (along each axis)
    #    - time
    f = h5py.File(os.path.join(dir_name,block_data[0][1]),'r')

    #print(f.attrs.keys())
    left_edge = (f.attrs['lower']) # ordered as x,y,z
    right_edge = (f.attrs['upper']) # ordered as x,y,z
    
    block_group = f[block_data[0][0]]
    ghost = (block_group.attrs['enzo_GridStartIndex'])
    start_index = (block_group.attrs['enzo_GridStartIndex']) # x,y,z
    end_index = (block_group.attrs['enzo_GridEndIndex']) #x,y,z
    # end_index is the last index included in the active zone
    res = (end_index+1-start_index) # x,y,z
    # Note that the actual fields are stored as arrays which are 
    # ordered as z,y,x
    time = block_group.attrs["time"]
    
    field_list = []
    for data_field_name in block_group.keys():
        if data_field_name[:6] == "field_":
            field_list.append((data_field_name,data_field_name[6:]))
        # do not presently support particle fields
    f.close()
    
    return res, field_list, left_edge, right_edge, ghost, time

def load_block_fields_(block_array, out_data, res, dir_name, 
                       field_list):
    # all 3 arguments are expected to be listed as (z,y,x)
    # face centered fields will not be properly loaded
    for k in range(block_array.shape[0]):
        for j in range(block_array.shape[1]):
            for i in range(block_array.shape[2]):
                
                block_name, block_file = block_array[k,j,i]
                
                data_idx = (slice(k*res[0],(k+1)*res[0]),
                            slice(j*res[1],(j+1)*res[1]),
                            slice(i*res[2],(i+1)*res[2]))
                
                f = h5py.File(os.path.join(dir_name,block_file),'r')
                block_group = f[block_name]
                block_idx = _get_slices(block_group.attrs)
                for save_name, out_name in field_list:
                    out_data[out_name][data_idx] = \
                        block_group[save_name][block_idx]
                f.close()

def load_full_grid(fname):
    """
    Used to load in the full static enzo-e grid at once (not AMR).
    
    For a full large simulation run, use yt for analysis. The main purpose 
    of this is for regression testing. It is motivated
        - The build system uses python 2. and the enzo_p yt frontend 
          encounters bugs for python 2.7
        - yt also doesn't support loading root meshes of [1,1,1]

    Outer Ghost zones get excluded

    Parameter
    ---------
    fname: string
        The path to the appropriate blocklist file

    Returns
    -------
    data_dict: dict
        Contains the arrays for each field. Each array is ordered as 
        [z][y][x]
    left_edge: np.ndarray, shape =(3,)
        An array of values (ordered as x,y,z) that indicates the 
        position of the upper left corner
    right_edge: np.ndarray, shape =(3,)
        An array of values (ordered as x,y,z) that indicates the 
        position of the lower right corner
    time: float
        Indicates the time of the simulation
    """
    
    assert os.path.isfile(fname)
    dir_name = os.path.dirname(fname)
    with open(fname, "r") as f:
        lines = f.readlines()
        # each entry in block_data is a tuple (block, block_file)
        block_data = [tuple(line.strip().split()) for line in lines]

    # root_shape is listed as (x,y,z)
    # block_array is shaped as (z,y,x)
    root_shape, block_array = _prepare_array_of_blocks(block_data)

    # all of the following multidimensional values are listed as (x,y,z)
    res, field_list, left_edge, right_edge, ghost, time = \
        get_details_(dir_name, block_data)

    # grid_shape is listed as (x,y,z)
    grid_shape = np.array(root_shape) * res

    data_dict = {}
    for save_name, out_name in field_list:
        data_dict[out_name] = np.zeros(shape = grid_shape[::-1],
                                       dtype = np.float64)
    load_block_fields_(block_array, data_dict, res[::-1], dir_name,
                       field_list)
    return data_dict, left_edge, right_edge, time[0]

def load_enzoe(fname, periodicity = (True, True, True)):
    # basically attempts to read enzo-e data into yt as a generic dataset
    
    
    # we are going to see what happens if we just tell it Cartesian and pass 
    # in x,y,z
    data_dict, left_edge, right_edge, time = load_full_grid(fname)
    
    known_fields = [('density', 'code_mass/code_length**3'),# "g/cm**3"),
                    ('velocity_x', 'code_velocity'),# "cm/s"),
                    ('velocity_y', 'code_velocity'),#, "cm/s"),
                    ('velocity_z', 'code_velocity'), #"cm/s"),
                    ('pressure', 'code_mass/code_length/code_time**2'),#"dyne/cm**2")
                    ("bfield_x", "code_magnetic"), #'gauss'
                    ("bfield_y", "code_magnetic"), #'gauss'
                    ("bfield_z", "code_magnetic"), #'gauss'
                   ]
    new_dict = {}
    shape = None
    for name, unit in known_fields:
        if name in data_dict:
            new_dict[name] = (data_dict[name], unit)
            shape = data_dict[name].shape
    for key in data_dict.keys():
        if key not in new_dict:
            new_dict[key] = data_dict
            shape = data_dict[name].shape

    bbox = np.array(list(zip(left_edge,right_edge)))
    ds = yt.load_uniform_grid(data=new_dict,domain_dimensions=shape,
                              bbox=bbox,sim_time = time, geometry = 'cartesian')
    return ds
