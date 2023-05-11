// See LICENSE_CELLO file for license and copyright information

/// @file     disk_utils.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2022-12-20
/// @brief    [\ref Disk] Declaration of utility functions used for debugging

// it may be better to place these functions in the array component
#include "_array.hpp"
#include <utility> // std::pair
#include <set>     // std::set

namespace disk_utils{


namespace detail{ // helper functions not meant to be used elsewhere

/// helper function that actually writes ``arr`` to an open HDF5 file
///
/// This assumes that ``file`` has been opened and that the appropriate group
/// has been created (and has already been opened)
template<typename T>
void dump_array_to_hdf5_(FileHdf5& file, const std::string &arr_name,
                         CelloArray<T,3> arr)
{
  const int type = cello::get_type_enum<T>();

  // for simplicity, we are going to take some steps to ensure that array
  // dimension and array size are identical (this is inefficient, but okay
  // since this is all for debugging purposes)
  int nx, ny, nz;
  nz = arr.shape(0);
  ny = arr.shape(1);
  nx = arr.shape(2);

  if (ny == 1){
    // I'm not totally sure why we run into problems when this happens
    ERROR("disk_utils::dump_array_to_hdf5_",
          "Something goes wrong when arr.shape(1) is 1.");
  }

  // allocate a contiguous array
  CelloArray<T,3> arr_copy(nz,ny,nx);
  // copy data into arr_copy
  arr.copy_to(arr_copy);

  file.mem_create(nx,ny,nz, nx,ny,nz, 0,0,0);

  file.data_create(arr_name.c_str(), type, nz,ny,nx,1, nz,ny,nx,1);
  file.data_write(arr_copy.data());
  file.data_close();

  file.mem_close();
}

} // namespace detail

inline std::string default_dump_group_name() { return "/data-dump"; }

/// create an HDF5 file called ``fname`` and then save copies of arrays to it
///
/// the array is saved in a group called "data-dump"
///
/// @param[in] fname is the file name where data will be saved
/// @param[in] pairs is a vector of (array-name, array) pairs that will be
///     saved to disk. No array-name should be duplicated
template<typename T>
void dump_arrays_to_hdf5(const std::string& fname,
                         const std::vector<std::pair<std::string,
                                           CelloArray<T,3>>> &pairs){
  // create the output file
  FileHdf5 file("./", fname);
  file.file_create();

  // create a file group to hold the data
  std::string group_name = default_dump_group_name();
  file.group_chdir(group_name);
  file.group_create();

  // initialze set to ensure names aren't duplicated
  std::set<std::string> name_set;

  for (const auto &pair: pairs){
    const std::string &arr_name = pair.first;
    CelloArray<T,3> arr = pair.second;

    if (name_set.find(arr_name) != name_set.end()){
      ERROR1("disk_utils::dump_arrays_to_hdf5",
             "More than 1 array was specified with the name %s",
             arr_name.c_str());
    } else {
      // insert name into name_set so we can ensure there aren't subsequent
      // name conflicts
      name_set.insert(arr_name);
    }

    // actually dump the array to the file
    detail::dump_array_to_hdf5_(file, arr_name, arr);
  }

  file.group_close();
  file.file_close();
}

/// create an HDF5 file called ``fname`` and then save a copy of ``arr`` to it
///
/// the array is saved in a group called "data-dump" with the name "data"
///
/// @param[in] fname is the file name where data will be saved
/// @param[in] arr is the array that is saved to disk
template<typename T>
void dump_array_to_hdf5(const std::string& fname, CelloArray<T,3> arr){
  std::vector<std::pair<std::string, CelloArray<T,3>>> tmp_v = {{"data", arr}};
  dump_arrays_to_hdf5(fname, tmp_v);
}

} // namespace disk_utils
