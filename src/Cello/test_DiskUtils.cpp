// See LICENSE_CELLO file for license and copyright information

/// @file     test_FileHdf5.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Fri Jan 20 2023
/// @brief    Program implementing unit tests for functions in the disk_utils
///           namespace

#include "main.hpp" 
#include "test.hpp"
#include "array.hpp"
#include "disk.hpp"

#include "test_ArrayTestTools.hpp" // range_3Darray, assert_allequal3D


namespace{ // anonymous namespace defines things that are only used locally

template<typename T>
std::string type_enum_name()
{ return cello::type_name[cello::get_type_enum<T>()]; }

/// read CelloArray from an HDF5 file
template<typename T>
CelloArray<T,3> read_data_(const std::string &fname,
                           const std::string &arr_name)
{
  FileHdf5 file("./",fname);
  file.file_open();

  int type = type_unknown;
  int mx, my, mz;
  file.data_open (disk_utils::default_dump_group_name() + "/" + arr_name,
                  &type, &mz, &my, &mx); // I don't know why we are reversing
                                         // order here, but it seems to work

  

  if (type != cello::get_type_enum<T>()){
    ERROR("read_data_", "There was a type mismatch!");
  }

  // allocate output
  CelloArray<T,3> out(mz, my, mx);

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
void test_dump_array_to_hdf5(){

  const std::string type_name = type_enum_name<T>();
  std::string func_name = "dump_array_to_hdf5[" + type_name + "]";
  unit_func(func_name.c_str());

  CkPrintf("Testing: %s\n", func_name.c_str());

  const std::string fname = "single-array-" + type_name + ".h5";

  // initialize an array with unique values
  CelloArray<T,3> arr = range_3Darray(4, 3, 2, (T)0, (T)1);

  // dump the array to disk
  disk_utils::dump_array_to_hdf5(fname, arr);

  // now load data from the file
  CelloArray<T,3> loaded_arr = read_data_<T>(fname, "data");

  // now confirm that the loaded data is identical to arr
  assert_allequal3D(arr, loaded_arr);

}

//----------------------------------------------------------------------

template<typename T>
void test_dump_arrays_to_hdf5(){

  const std::string type_name = type_enum_name<T>();
  std::string func_name = "dump_arrays_to_hdf5[" + type_name + "]";
  unit_func(func_name.c_str());

  CkPrintf("Testing: %s\n", func_name.c_str());

  const std::string fname = "multi-array-" + type_name + ".h5";

  // initialize an array of pairs with unique values
  std::vector<std::pair<std::string, CelloArray<T,3>>> orig_vec;

  orig_vec.push_back( { "array 1", range_3Darray(4, 3, 2, (T)0, (T)1) } );
  orig_vec.push_back( { "array 2", range_3Darray(2, 3, 4, (T)0, (T)-1) } );
  orig_vec.push_back( { "array 3", range_3Darray(4, 2, 1, (T)30, (T)2) } );
  orig_vec.push_back( { "array 4", range_3Darray(1, 2, 4, (T)30, (T)2) } );

  // dump the array to disk
  disk_utils::dump_arrays_to_hdf5(fname, orig_vec);

  // now load data from the file
  for (const auto &pair: orig_vec){
    CelloArray<T,3> loaded_arr = read_data_<T>(fname, pair.first);

    // now confirm that the loaded data is identical to arr
    assert_allequal3D(pair.second, loaded_arr);
  }
}

//----------------------------------------------------------------------

PARALLEL_MAIN_BEGIN
{
  PARALLEL_INIT;

  unit_init(0,1);

  test_dump_array_to_hdf5<float>();
  test_dump_array_to_hdf5<double>();

  test_dump_arrays_to_hdf5<float>();
  test_dump_arrays_to_hdf5<double>();

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END
