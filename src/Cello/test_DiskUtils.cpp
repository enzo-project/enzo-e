// See LICENSE_CELLO file for license and copyright information

/// @file     test_FileHdf5.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Fri Jan 20 2023
/// @brief    Program implementing unit tests for functions in the disk_utils
///           namespace

#include "main.hpp" 
#include "test.hpp"
#include "view.hpp"
#include "disk.hpp"

#include "test_ViewTestTools.hpp" // range_3Darray, assert_allequal3D


namespace{ // anonymous namespace defines things that are only used locally

template<typename T>
std::string type_enum_name()
{ return cello::type_name[cello::get_type_enum<T>()]; }

/// read CelloView from an HDF5 file
template<typename T>
CelloView<T,3> read_data_(const std::string &fname,
                           const std::string &view_name)
{
  FileHdf5 file("./",fname);
  file.file_open();

  int type = type_unknown;
  int mx, my, mz;
  file.data_open (disk_utils::default_dump_group_name() + "/" + view_name,
                  &type, &mz, &my, &mx); // I don't know why we are reversing
                                         // order here, but it seems to work

  

  if (type != cello::get_type_enum<T>()){
    ERROR("read_data_", "There was a type mismatch!");
  }

  // allocate output
  CelloView<T,3> out(mz, my, mx);

  // read the data
  file.data_read(out.data());

  // cleanup! (unclear how necessary this is)
  file.data_close();
  file.file_close();

  return out;
}

} // close namespace

//----------------------------------------------------------------------

template<typename T>
void test_dump_view_to_hdf5(){

  const std::string type_name = type_enum_name<T>();
  std::string func_name = "dump_view_to_hdf5[" + type_name + "]";
  unit_func(func_name.c_str());

  CkPrintf("Testing: %s\n", func_name.c_str());

  const std::string fname = "single-view-" + type_name + ".h5";

  // initialize an view with unique values
  CelloView<T,3> view = range_3Darray(4, 3, 2, (T)0, (T)1);

  // dump the view to disk
  disk_utils::dump_view_to_hdf5(fname, view);

  // now load data from the file
  CelloView<T,3> loaded_view = read_data_<T>(fname, "data");

  // now confirm that the loaded data is identical to view
  assert_allequal3D(view, loaded_view);

}

//----------------------------------------------------------------------

template<typename T>
void test_dump_views_to_hdf5(){

  const std::string type_name = type_enum_name<T>();
  std::string func_name = "dump_views_to_hdf5[" + type_name + "]";
  unit_func(func_name.c_str());

  CkPrintf("Testing: %s\n", func_name.c_str());

  const std::string fname = "multi-view-" + type_name + ".h5";

  // initialize a vector of pairs with unique values
  std::vector<std::pair<std::string, CelloView<T,3>>> orig_vec;

  orig_vec.push_back( { "view 1", range_3Darray(4, 3, 2, (T)0, (T)1) } );
  orig_vec.push_back( { "view 2", range_3Darray(2, 3, 4, (T)0, (T)-1) } );
  orig_vec.push_back( { "view 3", range_3Darray(4, 2, 1, (T)30, (T)2) } );
  orig_vec.push_back( { "view 4", range_3Darray(1, 2, 4, (T)30, (T)2) } );

  // dump the data to disk
  disk_utils::dump_views_to_hdf5(fname, orig_vec);

  // now load data from the file
  for (const auto &pair: orig_vec){
    CelloView<T,3> loaded_view = read_data_<T>(fname, pair.first);

    // now confirm that the loaded data is identical to view
    assert_allequal3D(pair.second, loaded_view);
  }
}

//----------------------------------------------------------------------

PARALLEL_MAIN_BEGIN
{
  PARALLEL_INIT;

  unit_init(0,1);

  test_dump_view_to_hdf5<float>();
  test_dump_view_to_hdf5<double>();

  test_dump_views_to_hdf5<float>();
  test_dump_views_to_hdf5<double>();

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END
