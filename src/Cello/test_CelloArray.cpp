// See LICENSE_CELLO file for license and copyright information

/// @file     test_CelloArray.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2020-06-09
/// @brief    Test program for the CelloArray class template
///
/// The main testing strategy here is to test some of the core functionallity
/// of CelloArray with:
///    - SimpleInitializationTests
///    - SimpleAssignmentTests
///    - SizeShapeTests
/// And then use more sophisticated machinery that relies on this
/// functionallity to test the remaining features.

#include <cstdlib>
#include <initializer_list>
#include <vector>

#include "main.hpp"
#include "test.hpp"
#include "array.hpp"

//----------------------------------------------------------------------

template<typename T, std::size_t D>
class MemManagedArrayBuilder{
  // This template class is used to build and a manage the lifetime of a
  // CelloArray that manages its own memory
  //
  // In the future it might be a good idea to have this and
  // PtrWrapArrayBuilder have a shared base class
public:
  template<typename... Args>
  MemManagedArrayBuilder(Args... args) : arr_(args...) {}
  CelloArray<T,D>* get_arr() { return &arr_; }
  T* get_wrapped_ptr() { return NULL; }
  static std::string name() { return "MemManagedArrayBuilder"; }

private:
  CelloArray<T,D> arr_;
};

//----------------------------------------------------------------------

template<typename T, std::size_t D>
class PtrWrapArrayBuilder{
  // This template class is used to build and a manage the lifetime of a
  // CelloArray wraps an existing pointer. This template class also manages the
  // lifetime of the underlying pointer.
  //
  // In the future it might be a good idea to have this and
  // MemManagedArrayBuilder have a shared base class
public:

  template<typename... Args>
  PtrWrapArrayBuilder(Args... args)
    : arr_ptr_(nullptr),
      ptr_(nullptr)
  {
    // we are sort of cheating here to compute the required size
    CelloArray<T,D> temp(args...);
    ptr_ = new T [temp.size()]{};
    arr_ptr_ = new CelloArray<T,D>(ptr_,args...);
  }

  ~PtrWrapArrayBuilder() {
    delete arr_ptr_;
    delete[] ptr_;
  }

  CelloArray<T,D>* get_arr() {return arr_ptr_;}
  T* get_wrapped_ptr() { return ptr_; }
  static std::string name() { return "PtrWrapArrayBuilder"; }

private:
  CelloArray<T,D>* arr_ptr_;
  T* ptr_;
};

//----------------------------------------------------------------------

void pointer_compare_(double vals[], std::vector<double> ref,
                      std::string func_name, const char* file, int line){
  bool all_match = true;
  for(std::size_t i = 0; i < ref.size(); i++){
    bool match = (ref[i] == vals[i]);
    if (!match){
      if (all_match){
        CkPrintf("\nUnequal Pointer Element Error in %s:\n",
                 func_name.c_str());
      }
      CkPrintf("Expected %e at index %lu. Got: %e\n", ref[i], i, vals[i]);
      all_match = false;
    }
  }
  Unit::instance()->assertion(all_match, file, line, true);
}

//----------------------------------------------------------------------

void compare_against_arr_(CelloArray<double, 2> &arr2d,
                          std::vector<double> ref,
                          std::string func_name,
                          const char* file, int line){

  int my = arr2d.shape(0);
  int mx = arr2d.shape(1);

  bool all_match = true;

  for(int iy=0; iy<my; iy++){
    for(int ix=0; ix<mx; ix++){
      int index = ix + mx*iy;
      bool match = (arr2d(iy,ix) == ref[index]);
      if (!match){
        if (all_match){
          CkPrintf("\nUnequal Array Element Error in %s:\n",
                   func_name.c_str());
        }
        all_match = false;
        CkPrintf("Expected %e at (%d, %d). Got: %e\n",
                 ref[index], iy, ix, arr2d(iy,ix));
      }
    }
  }
  Unit::instance()->assertion(all_match, file, line, true);
}

//----------------------------------------------------------------------

template<typename Builder>
void compare_builder_arr_(Builder &builder, std::vector<double> ref,
                          std::string func_name,
                          const char* file, int line){
  compare_against_arr_(*builder.get_arr(), ref, func_name, file, line);
  if (builder.get_wrapped_ptr() != nullptr){
    pointer_compare_(builder.get_wrapped_ptr(), ref, func_name, file, line);
  }
}

//----------------------------------------------------------------------

#define check_pointer_vals(VALS, REF, FUNC_NAME)                     \
  pointer_compare_(VALS, REF, FUNC_NAME, __FILE__, __LINE__);

#define check_arr_vals(VALS, REF, FUNC_NAME)                         \
  compare_against_arr_(VALS, REF, FUNC_NAME, __FILE__, __LINE__);

#define check_builder_arr(BUILDER, REF, FUNC_NAME)                   \
  compare_builder_arr_(BUILDER, REF, FUNC_NAME, __FILE__, __LINE__);

//----------------------------------------------------------------------

class SimpleInitializationTests{
  // This class holds tests used to check the simple initialization of
  // CelloArray. Unfortunately, it is not really possible to completely
  // disentangle the following set of features and test them completely
  // independently of each other:
  // 1. direct initialization of an instance that wraps an existing pointer
  // 2. direct initialization of an instance which manages its own data
  // 3. retrieval operations
public:

  void test_init_simple_managed_memory(){

    // all values are assumed to be zero-initialized
    std::string template_message_str =
      ("There is an issue in either the constructor for instance of %dD "
       "CelloArrays that manage own memory, or in CelloArray::operator()");

    const char *template_message = template_message_str.c_str();
    
    CelloArray<double, 1> arr1D(8);
    for (int i = 0; i<8; i++){
      ASSERT1("IrreducibleTests::test_managed_memory", template_message,
              1, arr1D(i) == 0);
    }


    CelloArray<double, 2> arr2D(2,4);
    for (int i = 0; i<2; i++){
      for (int j = 0; j<4; j++){
        ASSERT1("IrreducibleTests::test_managed_memory", template_message,
                2, arr2D(i,j) == 0);
      }
    }

    CelloArray<double, 3> arr3D(2,3,4);
    for (int i = 0; i<2; i++){
      for (int j = 0; j<3; j++){
        for (int k = 0; k<4; k++){
          ASSERT1("IrreducibleTests::test_managed_memory", template_message,
                  3, arr3D(i,j,k) == 0);
        }
      }
    }

    CelloArray<double, 4> arr4D(2,3,4,5);
    for (int i = 0; i<2; i++){
      for (int j = 0; j<3; j++){
        for (int k = 0; k<4; k++){
          for (int l = 0; l<5; l++){
            ASSERT1("IrreducibleTests::test_managed_memory", template_message,
                    4, arr4D(i,j,k,l) == 0);
          }
        }
      }
    }
  }

  void test_init_simple_wrapped_pointer(){
    // all values are assumed to be zero-initialized
    std::string template_message_str =
      ("There is an issue in either the constructor for instance of %dD "
       "CelloArrays that wraps a pointer, or in CelloArray::operator()");

    const char *template_message = template_message_str.c_str();

    double zero_values[120] = { };
    double* ptr = &(zero_values[0]);
    
    CelloArray<double, 1> arr1D(ptr, 8);
    for (int i = 0; i<8; i++){
      ASSERT1("IrreducibleTests::test_managed_memory", template_message,
              1, arr1D(i) == 0);
    }


    CelloArray<double, 2> arr2D(ptr,2,4);
    for (int i = 0; i<2; i++){
      for (int j = 0; j<4; j++){
        ASSERT1("IrreducibleTests::test_managed_memory", template_message,
                2, arr2D(i,j) == 0);
      }
    }

    CelloArray<double, 3> arr3D(ptr,2,3,4);
    for (int i = 0; i<2; i++){
      for (int j = 0; j<3; j++){
        for (int k = 0; k<4; k++){
          ASSERT1("IrreducibleTests::test_managed_memory", template_message,
                  3, arr3D(i,j,k) == 0);
        }
      }
    }


    CelloArray<double, 4> arr4D(ptr,2,3,4,5);
    for (int i = 0; i<2; i++){
      for (int j = 0; j<3; j++){
        for (int k = 0; k<4; k++){
          for (int l = 0; l<5; l++){
            ASSERT1("IrreducibleTests::test_managed_memory", template_message,
                    4, arr4D(i,j,k,l) == 0);
          }
        }
      }
    }
  }

  void test_init_nonzero_wrapped_pointer(){

    std::string template_message_str =
      ("There is an issue in either the constructor for instance of %dD "
       "CelloArrays that wraps a pointer, or in CelloArray::operator()");
    const char *template_message = template_message_str.c_str();

    double arr_vals[24] = { 0.,  1.,  2.,  3.,  4.,  5.,
                            6.,  7.,  8.,  9., 10., 11.,
                           12., 13., 14., 15., 16., 17.,
                           18., 19., 20., 21., 22., 23.};
    double* ptr = &(arr_vals[0]);

    CelloArray<double, 1> arr1D(ptr, 13);
    for (int i = 0; i<13; i++){
      ASSERT1("IrreducibleTests::test_managed_memory", template_message,
              1, arr1D(i) == arr_vals[i]);
    }


    CelloArray<double, 2> arr2D(ptr,3,5);
    for (int i = 0; i<3; i++){
      for (int j = 0; j<5; j++){
        ASSERT1("IrreducibleTests::test_managed_memory", template_message,
                2, arr2D(i,j) == arr_vals[j+i*5]);
      }
    }

    CelloArray<double, 3> arr3D(ptr,2,3,4);
    for (int i = 0; i<2; i++){
      for (int j = 0; j<3; j++){
        for (int k = 0; k<4; k++){
          ASSERT1("IrreducibleTests::test_managed_memory", template_message,
                  3, arr3D(i,j,k) == arr_vals[k+4*(j+3*i)]);
        }
      }
    }

  }

  void run_tests(){
    test_init_simple_managed_memory();
    // it is not possible to initialize a non-zero pointer
    test_init_simple_wrapped_pointer();
    test_init_nonzero_wrapped_pointer();
  }
};

//----------------------------------------------------------------------

class SimpleElementAssignmentTests{
  // Used to test whether individual element assignments actually work
private:
  
  std::string get_assignment_template_msg_(bool wrapped_pointer){
    std::string template_message_str;
    if (wrapped_pointer){
      template_message_str =
        ("There is an issue with operator() for instances of %dD CelloArrays "
         "that own their own memory");
    } else {
      template_message_str =
        ("There is an issue with operator() for instances of %dD CelloArrays "
         "that wrap existing pointers");
    }
    return template_message_str;
  }

  void test_assignment_1D_(unsigned int seed, double *wrapped_pointer,
                           CelloArray<double, 1> &arr1D, int length)
  {
    std::string template_message_str =
      get_assignment_template_msg_(wrapped_pointer != NULL);
    const char *template_message = template_message_str.c_str();

    int npasses = 2;
    if (wrapped_pointer != NULL) { npasses++;}

    for (int t = 0; t<npasses; t++){
      std::srand(seed);
      for (int i = 0; i<8; i++){

        double val = ((double)std::rand()/RAND_MAX);
        if (t == 0){
          arr1D(i) = val;
        } else if (t == 1){
          ASSERT1("IrreducibleTests::test_assignment_managed_memory",
                  template_message, 1, arr1D(i) == val);
        } else {
          ASSERT1("IrreducibleTests::test_assignment_managed_memory",
                  "The wrapped pointer of a 1D CelloArray is not updated",
                  1, wrapped_pointer[i] == val);
        }
      }
    }
  }

  void test_assignment_3D_(unsigned int seed, double *wrapped_pointer,
                           CelloArray<double, 3> &arr3D, int mz, int my,
                           int mx)
  {
    std::string template_message_str =
      get_assignment_template_msg_(wrapped_pointer != NULL);
    const char *template_message = template_message_str.c_str();

    int npasses = 2;
    if (wrapped_pointer != NULL) { npasses++;}

    for (int t = 0; t<npasses; t++){
      std::srand(seed);
      for (int iz = 0; iz<mz; iz++){
        for (int iy = 0; iy<my; iy++){
          for (int ix = 0; ix<mx; ix++){

            double val = ((double)std::rand()/RAND_MAX);
            if (t == 0){
              arr3D(iz,iy,ix) = val;
            } else if (t == 1){
              ASSERT1("IrreducibleTests::test_assignment_3D_",
                      template_message, 3, arr3D(iz,iy,ix) == val);
            } else {
              ASSERT("IrreducibleTests::test_assignment_3D_",
                     "The wrapped pointer of a 3D CelloArray has not been "
                     "appropriately updated",
                     wrapped_pointer[ix+mx*(iy+my*iz)] == val);
            }
          }
        }
      }
    }
  }

  void test_assignment_4D_(unsigned int seed, double *wrapped_pointer,
                           CelloArray<double, 4> &arr4D, int mz, int my,
                           int mx, int mw)
  {
    std::string template_message_str =
      get_assignment_template_msg_(wrapped_pointer != NULL);
    const char *template_message = template_message_str.c_str();

    int npasses = 2;
    if (wrapped_pointer != NULL) { npasses++;}

    for (int t = 0; t<npasses; t++){
      std::srand(seed);
      for (int iz = 0; iz<mz; iz++){
        for (int iy = 0; iy<my; iy++){
          for (int ix = 0; ix<mx; ix++){
            for (int iw = 0; iw<mw; iw++){

              double val = ((double)std::rand()/RAND_MAX);
              if (t == 0){
                arr4D(iz,iy,ix,iw) = val;
              } else if (t == 1){
                ASSERT1("IrreducibleTests::test_assignment_4D_",
                        template_message, 4, arr4D(iz,iy,ix,iw) == val);
              } else {
                ASSERT("IrreducibleTests::test_assignment_4D_",
                       "The wrapped pointer of a 4D CelloArray has not been "
                       "appropriately updated",
                       wrapped_pointer[iw+mw*(ix+mx*(iy+my*iz))] == val);
              }
            }
          }
        }
      }
    }
  }

  void test_assignment_managed_memory(){

    unsigned int seed = 5342342;
    CelloArray<double, 1> arr1D(8);
    test_assignment_1D_(seed, NULL, arr1D, 8);

    CelloArray<double, 3> arr3D(2,3,4);
    test_assignment_3D_(seed, NULL, arr3D, 2,3,4);

    CelloArray<double, 4> arr4D(2,3,4,5);
    test_assignment_4D_(seed, NULL, arr4D, 2,3,4,5);
  }

  void test_assignment_wrapped_ptr(){

    unsigned int seed = 5342342;
    double *ptr;

    double temp_arr1[120] = {};
    ptr = &(temp_arr1[0]);
    CelloArray<double, 1> arr1D(ptr,8);
    test_assignment_1D_(seed, ptr, arr1D, 8);

    double temp_arr3[120] = {};
    ptr = &(temp_arr3[0]);
    CelloArray<double, 3> arr3D(ptr,2,3,4);
    test_assignment_3D_(seed, ptr, arr3D, 2,3,4);

    double temp_arr4[120] = {};
    ptr = &(temp_arr4[0]);
    CelloArray<double, 4> arr4D(ptr,2,3,4,5);
    test_assignment_4D_(seed, ptr, arr4D, 2,3,4,5);
  }

public:
  void run_tests(){
    test_assignment_managed_memory();
    test_assignment_wrapped_ptr();
  }
};

//----------------------------------------------------------------------

class SizeShapeTests{
  // Used to test the size and shape methods of CelloArray
private:

  template<std::size_t D>
  void check_expected_size_shape_(CelloArray<double, D> &arr,
                                  std::initializer_list<int> expect_shape){
    static_assert(D>0, "D must always be positive.");
    ASSERT("SizeShapeTests::check_expected_size_shape_",
           "the length of the expect_shape initializer_list is wrong",
           expect_shape.size() == D);

    std::size_t axis_ind = 0;
    std::size_t cur_size = 1;
    for (const int axis_length : expect_shape){
      ASSERT("SizeShapeTests::check_expected_size_shape_",
             ("the expect_shape initializer_list should only include positive "
              "values"), axis_length > 0);
      ASSERT3("SizeShapeTests::check_expected_size_shape_",
              "the length along axis %d was %d. It should be %d.",
              (int)axis_ind, arr.shape(axis_ind), axis_length,
              ((int)arr.shape(axis_ind)) == ((int)axis_length));
      cur_size *= (std::size_t)axis_length;
      axis_ind++;
    }

    ASSERT("SizeShapeTests::check_expected_size_shape_",
           "did not get the expected result for size",
           (std::size_t)arr.size() == (std::size_t)cur_size);
  }

  void test_size_shape_1D_(double* ptr, int length){
    if (ptr != NULL){
      CelloArray<double, 1> arr(ptr,length);
      check_expected_size_shape_(arr,{length});
    } else {
      CelloArray<double, 1> arr(length);
      check_expected_size_shape_(arr,{length});
    }
  }

  void test_size_shape_2D_(double* ptr, int my, int mx){
    if (ptr != NULL){
      CelloArray<double, 2> arr(ptr,my,mx);
      check_expected_size_shape_(arr,{my,mx});
    } else {
      CelloArray<double, 2> arr(my,mx);
      check_expected_size_shape_(arr,{my,mx});
    }
  }

  void test_size_shape_3D_(double* ptr, int mz, int my, int mx){
    if (ptr != NULL){
      CelloArray<double, 3> arr(ptr,mz,my,mx);
      check_expected_size_shape_(arr,{mz,my,mx});
    } else {
      CelloArray<double, 3> arr(mz,my,mx);
      check_expected_size_shape_(arr,{mz,my,mx});
    }
  }

  void test_size_shape_4D_(double* ptr, int mz, int my, int mx, int mw){
    if (ptr != NULL){
      CelloArray<double, 4> arr(ptr,mz,my,mx,mw);
      check_expected_size_shape_(arr,{mz,my,mx,mw});
    } else {
      CelloArray<double, 4> arr(mz,my,mx,mw);
      check_expected_size_shape_(arr,{mz,my,mx,mw});
    }
  }

  void test_size_shape_5D_(double* ptr, int mz, int my, int mx, int mw, int mv){
    if (ptr != NULL){
      CelloArray<double, 5> arr(ptr,mz,my,mx,mw,mv);
      check_expected_size_shape_(arr,{mz,my,mx,mw,mv});
    } else {
      CelloArray<double, 5> arr(mz,my,mx,mw,mv);
      check_expected_size_shape_(arr,{mz,my,mx,mw,mv});
    }
  }

  void test_size_shape_managed_(bool wrapped){
    double carray[750] = { };
    double *ptr = (wrapped) ? &(carray[0]) : NULL;
    test_size_shape_1D_(ptr, 1);
    test_size_shape_1D_(ptr, 24);

    test_size_shape_2D_(ptr,  1, 1);
    test_size_shape_2D_(ptr, 24, 1);
    test_size_shape_2D_(ptr,  1, 24);

    test_size_shape_3D_(ptr,  1,  1,  1);
    test_size_shape_3D_(ptr, 24,  1,  1);
    test_size_shape_3D_(ptr,  1, 58,  1);
    test_size_shape_3D_(ptr,  1,  1,  7);
    test_size_shape_3D_(ptr, 83,  5,  1);
    test_size_shape_3D_(ptr, 26,  1, 15);
    test_size_shape_3D_(ptr, 1,  42, 11);
    test_size_shape_3D_(ptr, 2,   3,  4);

    test_size_shape_4D_(ptr, 1,   1,  1,  1);
    test_size_shape_4D_(ptr, 2,   3,  4,  5);

    test_size_shape_5D_(ptr, 1,   1,  1,  1, 6);
    test_size_shape_5D_(ptr, 2,   3,  4,  5, 6);
  }

public:
  void run_tests(){
    test_size_shape_managed_(false);
    test_size_shape_managed_(true);
  }
};

//----------------------------------------------------------------------

class VariableAssignmentTests{
  // these are tests that check that check the assignments of arrays to
  // variables.

  template<template<typename, std::size_t> class Builder>
  void test_assignment_(){
    // vector used to hold expected results (the contents of this vector will
    // be updated throughout this test
    std::vector<double> expected;

    // construct vectors that hold pointers to shared copies of the tested
    // array to simplify the code required to makes sure that all of the
    // instances reflect the fact that they are shared copies
    std::vector<CelloArray<double, 2>*> array_ptr_vec;

    // vectors holding pointers to c-style arrays that at one time or another
    // were wrapped by a CelloArray
    std::vector<double*> wrapped_ptr_vec;
    // holds pointers to vectors containing the expected values that should be
    // held by each pointer that was wrapped at one time or another
    std::vector<std::vector<double>*> expected_wrapped_ptr_vals;

    // define a lambda function to actually carry out the checks on the
    // for the expected array values
    auto check_shared_copy_ = [&expected, &array_ptr_vec,
                               &wrapped_ptr_vec, &expected_wrapped_ptr_vals]()
      {
        for (CelloArray<double,2>* cur_arr_ptr : array_ptr_vec){
          check_arr_vals(*cur_arr_ptr, expected,
                         "VariableAssignmentTests::test_assignment_");
        }
        for (std::size_t i =0; i<wrapped_ptr_vec.size(); i++){
          check_pointer_vals(wrapped_ptr_vec[i],
                             *(expected_wrapped_ptr_vals[i]),
                             "VariableAssignmentTests::test_assignment_");
        }
      };

    // now let's set up the initial instance of CelloArray and initialize it
    Builder<double, 2> builder(2,3);
    CelloArray<double, 2> *arr_ptr = builder.get_arr();
    (*arr_ptr)(0,2) = 97.25;
    expected.assign({0, 0, 97.25,
                     0, 0,     0});

    // register it in array_ptr_vec and the wrapped array (if applicable) in
    // wrapped_ptr_vec
    array_ptr_vec.push_back(arr_ptr);
    if (builder.get_wrapped_ptr() != nullptr){
      wrapped_ptr_vec.push_back(builder.get_wrapped_ptr());
      expected_wrapped_ptr_vals.push_back(&expected);
    }
    // sanity check to confirm that all values are what we would expect
    check_shared_copy_();

    // Basically, this test consists of defining a collection of variables
    // holding CelloArrays that are each defined in a different way. The
    // variable is then assigned a shallow copy of the original array
    // Then the contents of the underlying values are modified and the contents
    // of all shallow copies are checked to make sure they are accurate

    // First, consider an initially default-initialized variable
    CelloArray<double, 2> second_var; // this is default constructed
    // assign the array pointed to by arr_ptr to second_var
    second_var = *arr_ptr;
    array_ptr_vec.push_back(&second_var);
    // confirm that the copy worked
    //print_array(second_var);
    check_shared_copy_();
    // check that they are all actually shallow copies
    second_var(0,1) = 4324.5;
    (*arr_ptr)(1,0) = -1;
    expected.assign({ 0, 4324.5, 97.25,
                     -1,     0.,    0.});
    // confirm that the copy worked
    check_shared_copy_();


    // Next, consider an initially defined array with the same shape
    CelloArray<double,2> third_var(2,3);
    third_var(1,2) = 42657.5;
    // now actually assign the value
    third_var = second_var;
    array_ptr_vec.push_back(&third_var);
    check_shared_copy_();

    third_var(1,2)  = -43.75;
    second_var(0,2) = 1500.;
    expected.assign({ 0, 4324.5,  1500.,
                     -1,    0.,  -43.75});
    // confirm that the copy worked
    check_shared_copy_();

    // Next, consider an initially defined array with a different shape
    CelloArray<double,2> fourth_var(100,53);
    fourth_var = third_var;
    array_ptr_vec.push_back(&fourth_var);
    check_shared_copy_();

    third_var(1,1)  = 12;
    expected.assign({ 0, 4324.5,  1500.,
                     -1,    12., -43.75});
    check_shared_copy_();

    // Next, consider an initially defined array that wraps a pointer
    double dummy_carray[6] = {0., 1., 2., 3., 4., 5.};
    double* dummy_wrapped_ptr = &(dummy_carray[0]);
    std::vector<double> expected_dummy_carray_vals(dummy_wrapped_ptr,
                                                   dummy_wrapped_ptr + 6);
    // sanity check
    CelloArray<double,2> fifth_var(dummy_wrapped_ptr, 2,3);
    array_ptr_vec.push_back(&fifth_var);
    wrapped_ptr_vec.push_back(dummy_wrapped_ptr);
    expected_wrapped_ptr_vals.push_back(&expected_dummy_carray_vals);
    fifth_var = *arr_ptr;
    // make sure that all values are expected!
    check_shared_copy_();

    // make sure that all shallow copies remain consistent after modification
    second_var(0,0)  = -546734;
    expected.assign({ -546734, 4324.5,  1500.,
                           -1,    12., -43.75});
    check_shared_copy_();
    // finally make sure that modifications to the wrapped vector doesn't
    // change values in each of the shallow copies of the original array
    dummy_carray[3] = 526;
    expected_dummy_carray_vals[3] = 526;
    check_shared_copy_();

    // Things that could still be checked:
    // 1. What happens when we deallocate an array (whether it's the root one
    //    or not)? [for the root copy case, when the array wraps an existing
    //    pointer, it would be interesting to ensure that the underlying
    //    pointer is modified] 
    // 2. What happens when we assign an array to another array (similar to
    //    deallocation questions)
    // 3. What happens when we assign an existing shallow copy equal to another
    //    shallow copy of the same underlying array?

  }

  template<template<typename, std::size_t> class Builder1,
           template<typename, std::size_t> class Builder2>
  void test_deepcopy_assignment_(){
    std::string func_name =
      ("VariableAssignmentTests::test_deepcopy_assignment_<" +
       Builder1<double, 2>::name() + "," +
       Builder2<double, 2>::name() + ">");
      
    Builder1<double, 2> builder_a(2,3);
    CelloArray<double, 2> *arr_ptr_a = builder_a.get_arr();
    (*arr_ptr_a)(0,2) = 97.25;

    Builder2<double, 2> builder_b(4,4);
    CelloArray<double, 2> *arr_ptr_b = builder_b.get_arr();
    (*arr_ptr_b)(1,1) = 15;

    *arr_ptr_a = arr_ptr_b->deepcopy();

    (*arr_ptr_a)(2,3) = -42;
    (*arr_ptr_b)(3,2) = 42;

    check_arr_vals(*arr_ptr_a, std::vector<double>({0,  0, 0,   0,
                                                    0, 15, 0,   0,
                                                    0,  0, 0, -42,
                                                    0,  0, 0,   0}),
                   func_name);
    if (builder_a.get_wrapped_ptr() != nullptr){
      check_pointer_vals(builder_a.get_wrapped_ptr(),
                         std::vector<double>({0,0,97.25, 0, 0,   0}),
                         func_name);
    }
    
    check_builder_arr(builder_b,
                      std::vector<double>({0,  0,  0, 0,
                                           0, 15,  0, 0,
                                           0,  0,  0, 0,
                                           0,  0, 42, 0}), func_name);
  }

public:
  void run_tests(){
    test_assignment_<MemManagedArrayBuilder>();
    test_assignment_<PtrWrapArrayBuilder>();
    test_deepcopy_assignment_<MemManagedArrayBuilder,   PtrWrapArrayBuilder>();
    test_deepcopy_assignment_<   PtrWrapArrayBuilder,MemManagedArrayBuilder>();
    test_deepcopy_assignment_<MemManagedArrayBuilder,MemManagedArrayBuilder>();
    test_deepcopy_assignment_<   PtrWrapArrayBuilder,   PtrWrapArrayBuilder>();
  }
};

//----------------------------------------------------------------------

class BulkAssignmentTest{
  // this is inspired by some tests that indicated that there were problems
  // with the bulk assignment machinery

public:

  template<template<typename, std::size_t> class Builder>
  void test_assign_from_scalar_(){
    std::string func_name = "BulkAssignmentTest::test_assign_from_scalar_";

    Builder<double, 2> builder(2,3);
    CelloArray<double, 2> *arr_ptr = builder.get_arr();
    (*arr_ptr)(0,2) = 97.25;

    (*arr_ptr).subarray() = 57.24;

    check_builder_arr(builder, std::vector<double>({57.24, 57.24, 57.24,
                                                    57.24, 57.24, 57.24}),
                      func_name);

    (*arr_ptr).subarray(CSlice(0,2), CSlice(0,-1)) = -3;
    check_builder_arr(builder, std::vector<double>({-3, -3, 57.24,
                                                    -3, -3, 57.24}),
                      func_name);

    CelloArray<double, 2> subarray = arr_ptr->subarray(CSlice(0,2),
                                                       CSlice(1,3));
    subarray.subarray() = 15;
    check_arr_vals(subarray,  std::vector<double>({15,15,15,15}), func_name);
    check_builder_arr(builder, std::vector<double>({-3, 15, 15,
                                                    -3, 15, 15}), func_name);
  }
  

  template<template<typename, std::size_t> class Builder1,
           template<typename, std::size_t> class Builder2>
  void test_assign_from_array_(){
    std::string func_name =
      ("BulkAssignmentTest::test_assign_from_array_<" +
       Builder1<double, 2>::name() + "," +
       Builder2<double, 2>::name() + ">");
    
    Builder1<double, 2> builder_1(2,3);
    CelloArray<double, 2> *arr_ptr_a = builder_1.get_arr();
    (*arr_ptr_a)(0,2) = 97.25;

    Builder2<double, 2> builder_2a(2,3);
    CelloArray<double, 2> *arr_ptr_2a = builder_2a.get_arr();
    Builder2<double, 2> builder_2b(2,3);
    CelloArray<double, 2> *arr_ptr_2b = builder_2b.get_arr();
    double val_2a = 0;
    double val_2b = -1;
    for (int iy = 0; iy<2; iy++){
      for (int ix = 0; ix <3; ix++){
        (*arr_ptr_2a)(iy,ix) = val_2a;
        (*arr_ptr_2b)(iy,ix) = val_2b;
        val_2a++;
        val_2b--;
      }
    }

    // sanity checks!
    check_builder_arr(builder_2a, std::vector<double>({0, 1, 2,
                                                       3, 4, 5}), func_name);
    check_builder_arr(builder_2b, std::vector<double>({-1, -2, -3,
                                                       -4, -5, -6}),
                      func_name);

    // now lets try an assignment of a full array
    (*arr_ptr_a).subarray() = *arr_ptr_2a;
    (*arr_ptr_a)(1,1)= -36;
    check_builder_arr(builder_1, std::vector<double>({0,   1, 2,
                                                      3, -36, 5}), func_name);
    // make sure array_2a is unaffected
    check_builder_arr(builder_2a, std::vector<double>({0, 1, 2,
                                                       3, 4, 5}), func_name);

    // now lets try another assignment
    (*arr_ptr_a).subarray() = *arr_ptr_2b;
    (*arr_ptr_a)(0,0)= 5;
    check_builder_arr(builder_1, std::vector<double>({ 5, -2, -3,
                                                      -4, -5, -6}), func_name);
    // make sure array_2a is unaffected
    check_builder_arr(builder_2b, std::vector<double>({-1, -2, -3,
                                                       -4, -5, -6}),
                      func_name);
  }


  template<template<typename, std::size_t> class Builder1,
           template<typename, std::size_t> class Builder2>
  void test_assign_subarrays_(){
    // These tests actually uncovered a longstanding bug
    std::string func_name =
      ("BulkAssignmentTest::test_assign_subarrays_<" +
       Builder1<double, 2>::name() + "," +
       Builder2<double, 2>::name() + ">");
    Builder1<double, 2> builder_1(4,4);
    CelloArray<double, 2> *arr_ptr_a = builder_1.get_arr();
    (*arr_ptr_a)(2,2) = 97.25;
    CelloArray<double, 2> subarray_1 = arr_ptr_a->subarray(CSlice(1,3),
                                                           CSlice(1,3));
    subarray_1(0,1) = 5;

    // sanity check:
    check_arr_vals(subarray_1, std::vector<double>({0,     5,
                                                    0, 97.25}), func_name);
    check_builder_arr(builder_1,  std::vector<double>({0, 0,     0, 0,
                                                       0, 0,     5, 0,
                                                       0, 0, 97.25, 0,
                                                       0, 0,     0, 0}),
                      func_name);

    Builder2<double, 2> builder_2(4,5);
    CelloArray<double, 2> *arr_ptr_2 = builder_2.get_arr();
    CelloArray<double, 2> subarray_2 = arr_ptr_2->subarray(CSlice(1,3),
                                                           CSlice(1,4));
    double val_2 = 0;
    for (int iy = 0; iy<4; iy++){
      for (int ix = 0; ix <5; ix++){
        (*arr_ptr_2)(iy,ix) = val_2;
        val_2++;
      }
    }

    // sanity check:
    check_builder_arr(builder_2,  std::vector<double>({ 0,  1,  2,  3,  4,
                                                        5,  6,  7,  8,  9,
                                                        10, 11, 12, 13, 14,
                                                        15, 16, 17, 18, 19}),
                      func_name);
    check_arr_vals(subarray_2, std::vector<double>({ 6,  7,  8,
                                                    11, 12, 13}), func_name);

    // NOW TO ACTUALLY PERFORM A TEST
    subarray_1.subarray() = subarray_2.subarray(CSlice(0,2),
                                                CSlice(0,2));

    check_arr_vals(subarray_1,  std::vector<double>({ 6,  7,
                                                     11, 12}), func_name);
    check_builder_arr(builder_1, std::vector<double>({0,  0,  0, 0,
                                                      0,  6,  7, 0,
                                                      0, 11, 12, 0,
                                                      0,  0,  0, 0}),
                      func_name);
    // extra sanity checks!
    subarray_1(0,0) = 97;
    check_arr_vals(subarray_1,  std::vector<double>({97,  7,
                                                     11, 12}), func_name);
    check_builder_arr(builder_1, std::vector<double>({0,  0,  0, 0,
                                                      0, 97,  7, 0,
                                                      0, 11, 12, 0,
                                                      0,  0,  0, 0}),
                      func_name);
    check_builder_arr(builder_2, std::vector<double>({ 0,  1,  2,  3,  4,
                                                       5,  6,  7,  8,  9,
                                                       10, 11, 12, 13, 14,
                                                       15, 16, 17, 18, 19}),
                     func_name);
    check_arr_vals(subarray_2,  std::vector<double>({ 6,  7,  8,
                                                     11, 12, 13}), func_name);

    // SECOND TEST:
    subarray_1.subarray(CSlice(0,2),
                        CSlice(0,2)) = subarray_2.subarray(CSlice(0,2),
                                                           CSlice(1,3));
    check_arr_vals(subarray_1, std::vector<double>({ 7,  8,
                                                    12, 13}), func_name);
    check_builder_arr(builder_1, std::vector<double>({0,  0,  0, 0,
                                                      0,  7,  8, 0,
                                                      0, 12, 13, 0,
                                                      0,  0,  0, 0}),
                      func_name);
    
  }

  void run_tests(){
    test_assign_from_scalar_<MemManagedArrayBuilder>();
    test_assign_from_scalar_<PtrWrapArrayBuilder>();

    test_assign_from_array_<MemManagedArrayBuilder,   PtrWrapArrayBuilder>();
    test_assign_from_array_<   PtrWrapArrayBuilder,MemManagedArrayBuilder>();
    test_assign_from_array_<MemManagedArrayBuilder,MemManagedArrayBuilder>();
    test_assign_from_array_<   PtrWrapArrayBuilder,   PtrWrapArrayBuilder>();

    test_assign_subarrays_<MemManagedArrayBuilder,MemManagedArrayBuilder>();
    test_assign_subarrays_<MemManagedArrayBuilder,   PtrWrapArrayBuilder>();
    test_assign_subarrays_<   PtrWrapArrayBuilder,MemManagedArrayBuilder>();
    test_assign_subarrays_<   PtrWrapArrayBuilder,   PtrWrapArrayBuilder>();
  }


};

//----------------------------------------------------------------------

class PassByValueTests{

private:

  // an array of 0 should be passed to this function
  void pass_by_val(CelloArray<double,2> arr,std::string func_name){
    // sanity check:
    check_arr_vals(arr, std::vector<double>({ 0, 0, 0,
                                              0, 0, 0}), func_name.c_str());
    arr(0,1) = 1;
    check_arr_vals(arr, std::vector<double>({ 0, 1, 0,
                                              0, 0, 0}), func_name.c_str());
  }
public:

  template<template<typename, std::size_t> class Builder>
  void test_pass_by_val_(){
    Builder<double, 2> builder(2,3);
    CelloArray<double, 2> *arr_ptr = builder.get_arr();

    pass_by_val(*arr_ptr, "PassByValueTests::test_pass_by_val");
    check_builder_arr(builder, std::vector<double>({ 0, 1, 0,
                                                     0, 0, 0}),
      "PassByValueTests::test_pass_by_val");
    
  }

  void run_tests(){
    test_pass_by_val_<MemManagedArrayBuilder>();
    test_pass_by_val_<PtrWrapArrayBuilder>();
  }

};


//----------------------------------------------------------------------

class IsAliasTests{

public:


  void test_is_alias_independent_ptr_wrapping_(){
    PtrWrapArrayBuilder<double, 2> builder(2,3);
    CelloArray<double, 2> *arr_ptr = builder.get_arr();

    double *ptr = builder.get_wrapped_ptr();

    CelloArray<double, 2> array_2(ptr + 3, 1, 3);

    ASSERT("IsAliasTests::test_is_alias_independent_ptr_wrapping_",
           "The arrays aren't aliases because they only partially overlap.",
           !arr_ptr->is_alias(array_2));

    CelloArray<double, 2> subarray = arr_ptr->subarray(CSlice(1,2),
                                                       CSlice(0,3));
    ASSERT("IsAliasTests::test_is_alias_independent_ptr_wrapping_",
           ("The arrays are aliases even though they were constructed "
            "independently from each other"),
           subarray.is_alias(array_2));

    CelloArray<double, 1> array_3(ptr + 3, 3);
    ASSERT("IsAliasTests::test_is_alias_independent_ptr_wrapping_",
           ("The arrays are not perfect aliases because they have different "
            "numbers of dimensions."), !subarray.is_alias(array_3));
  }
  
  template<template<typename, std::size_t> class Builder>
  void test_is_alias_(){
    Builder<double, 2> builder(2,3);
    CelloArray<double, 2> *arr_ptr = builder.get_arr();

    double data[6] = {0, 0, 0, 0, 0, 0};
    CelloArray<double, 2> array_2(data, 2, 3);

    ASSERT("IsAliasTests::test_is_alias_", "The arrays aren't aliases.",
           !arr_ptr->is_alias(array_2));

    CelloArray<double, 2> array_3(2, 3);

    ASSERT("IsAliasTests::test_is_alias_", "The arrays aren't aliases.",
           !arr_ptr->is_alias(array_3));

    CelloArray<double, 2> subarray1 = arr_ptr->subarray(CSlice(0,2),
                                                        CSlice(1,3));
    CelloArray<double, 2> subarray2 = arr_ptr->subarray(CSlice(0,2),
                                                        CSlice(0,2));
    ASSERT("IsAliasTests::test_is_alias_",
           "The arrays have partial but are not aliases.",
           !arr_ptr->is_alias(subarray1));
    ASSERT("IsAliasTests::test_is_alias_",
           "The arrays have partial but are not aliases.",
           !arr_ptr->is_alias(subarray2));

    CelloArray<double, 2> subarray1a = subarray1.subarray(CSlice(0,2),
                                                          CSlice(0,1));
    CelloArray<double, 2> subarray2a = subarray2.subarray(CSlice(0,2),
                                                          CSlice(1,2));
    ASSERT("IsAliasTests::test_is_alias_", "The arrays are aliases.",
           subarray1a.is_alias(subarray2a));
    ASSERT("IsAliasTests::test_is_alias_", "The arrays are aliases.",
           subarray1a.is_alias(arr_ptr->subarray(CSlice(0,2),
                                                 CSlice(1,2)))
           );
  }

  void run_tests(){
    test_is_alias_<MemManagedArrayBuilder>();
    test_is_alias_<PtrWrapArrayBuilder>();
    test_is_alias_independent_ptr_wrapping_();
  }

};

//----------------------------------------------------------------------

PARALLEL_MAIN_BEGIN
{
  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("CelloArray");
  
  SimpleInitializationTests init_tests;
  init_tests.run_tests();

  SimpleElementAssignmentTests element_assign_tests;
  element_assign_tests.run_tests();

  SizeShapeTests size_shape_tests;
  size_shape_tests.run_tests();

  VariableAssignmentTests var_assign_tests;
  var_assign_tests.run_tests();

  BulkAssignmentTest bulk_assignment_tests;
  bulk_assignment_tests.run_tests();

  PassByValueTests pass_by_val_tests;
  pass_by_val_tests.run_tests();

  IsAliasTests is_alias_tests;
  is_alias_tests.run_tests();

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END
