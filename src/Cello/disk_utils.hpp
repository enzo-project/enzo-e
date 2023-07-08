// See LICENSE_CELLO file for license and copyright information

/// @file     disk_utils.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2022-12-20
/// @brief    [\ref Disk] Declaration of utility functions used for debugging

// it may be better to place these functions in the view component
#include "_view.hpp"
#include <utility> // std::pair
#include <set>     // std::set

namespace disk_utils{


namespace detail{ // helper functions not meant to be used elsewhere

/// helper function that actually writes ``view`` to an open HDF5 file
///
/// This assumes that ``file`` has been opened and that the appropriate group
/// has been created (and has already been opened)
template<typename T>
void dump_view_to_hdf5_(FileHdf5& file, const std::string &view_name,
                        CelloView<T,3> view)
{
  const int type = cello::get_type_enum<T>();

  // for simplicity, we are going to take some steps to ensure that view
  // dimension and view size are identical (this is inefficient, but okay
  // since this is all for debugging purposes)
  int nx, ny, nz;
  nz = view.shape(0);
  ny = view.shape(1);
  nx = view.shape(2);

  if (ny == 1){
    // I'm not totally sure why we run into problems when this happens
    ERROR("disk_utils::dump_view_to_hdf5_",
          "Something goes wrong when view.shape(1) is 1.");
  }

  // allocate a contiguous view
  CelloView<T,3> view_copy(nz,ny,nx);
  // copy data into arr_copy
  view.copy_to(view_copy);

  file.mem_create(nx,ny,nz, nx,ny,nz, 0,0,0);

  file.data_create(view_name.c_str(), type, nz,ny,nx,1, nz,ny,nx,1);
  file.data_write(view_copy.data());
  file.data_close();

  file.mem_close();
}

} // namespace detail

inline std::string default_dump_group_name() { return "/data-dump"; }

/// create an HDF5 file called ``fname`` and then save copies of views to it
///
/// the view is saved in a group called "data-dump"
///
/// @param[in] fname is the file name where data will be saved
/// @param[in] pairs is a vector of (view-name, view) pairs that will be
///     saved to disk. No view-name should be duplicated
template<typename T>
void dump_views_to_hdf5(const std::string& fname,
                        const std::vector<std::pair<std::string,
                                                    CelloView<T,3>>> &pairs){
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
    const std::string &view_name = pair.first;
    CelloView<T,3> view = pair.second;

    if (name_set.find(view_name) != name_set.end()){
      ERROR1("disk_utils::dump_views_to_hdf5",
             "More than 1 view was specified with the name %s",
             view_name.c_str());
    } else {
      // insert name into name_set so we can ensure there aren't subsequent
      // name conflicts
      name_set.insert(view_name);
    }

    // actually dump the view to the file
    detail::dump_view_to_hdf5_(file, view_name, view);
  }

  file.group_close();
  file.file_close();
}

/// create an HDF5 file called ``fname`` and then save a copy of ``view`` to it
///
/// the view is saved in a group called "data-dump" with the name "data"
///
/// @param[in] fname is the file name where data will be saved
/// @param[in] view is the view that is saved to disk
template<typename T>
void dump_view_to_hdf5(const std::string& fname, CelloView<T,3> view){
  std::vector<std::pair<std::string, CelloView<T,3>>> tmp_v = {{"data", view}};
  dump_views_to_hdf5(fname, tmp_v);
}

} // namespace disk_utils
