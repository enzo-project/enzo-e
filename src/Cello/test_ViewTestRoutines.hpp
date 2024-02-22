// See LICENSE_CELLO file for license and copyright information

/// @file     test_ViewTestRoutines.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-07-03
/// @brief    Common routines used for testing CelloView and related machinery

#ifndef TEST_VIEW_TEST_ROUNTINES_HPP
#define TEST_VIEW_TEST_ROUNTINES_HPP

#include "test.hpp"
#include "view.hpp"
#include "error.hpp"

#include <array>
#include <vector>

// older test-code doesn't expect unit_assert to be used
#ifndef VIEW_TEST_ROUTINES_EMPLOY_UNIT_ASSERT
#define VIEW_TEST_ROUTINES_EMPLOY_UNIT_ASSERT 0
#endif

#if VIEW_TEST_ROUTINES_EMPLOY_UNIT_ASSERT != 0
#define VIEW_TEST_ASSERT(F,M,A)                                          \
  { if (!(A)) {                                                         \
      cello::message(stderr, "ERROR", __FILE__, __LINE__, F, M);         \
      unit_assert(false);                                                \
    } else { unit_assert(true); } }
#else
#define VIEW_TEST_ASSERT(F,M,A) ASSERT(F, M, A)
#endif

template<typename T>
struct ValRange{
  ValRange(T start, T step) : start(start), step(step) {}
  
  T start;
  T step;
};

//----------------------------------------------------------------------

template<typename T>
void assign_range_3Darray(const CelloView<T,3>& arr, ValRange<T> val_range) {
  int count = 0;

  const int mz = arr.shape(0);
  const int my = arr.shape(1);
  const int mx = arr.shape(2);

  for (int iz = 0; iz < mz; iz++){
    for (int iy = 0; iy < my; iy++){
      for (int ix = 0; ix < mx; ix++){
        arr(iz,iy,ix) = val_range.start + static_cast<T>(count)*val_range.step;
        count++;
      }
    }
  }

}

template<typename T>
void assign_range_3Darray(const CelloView<T,3>& arr, T start, T step) {
  assign_range_3Darray(arr, ValRange<T>(start, step));
}

//----------------------------------------------------------------------

template<typename T>
CelloView<T,3> range_3Darray(int mz, int my, int mx, ValRange<T> val_range) {
  CelloView<T,3> out(mz,my,mx);
  assign_range_3Darray(out, val_range);
  return out;
}

template<typename T>
CelloView<T,3> range_3Darray(int mz, int my, int mx, T start, T step) {
  return range_3Darray(mz, my, mx, ValRange<T>(start, step));
}

//----------------------------------------------------------------------

template<typename T>
void assert_allequal3D_(const CelloView<T,3>& a, const CelloView<T,3>& b,
                        const char* file, int line) {
  ASSERT6("assert_allequal3D",
          "The arrays don't have the same shape. Arrays a and b have shapes "
          "(%d,%d,%d) and (%d,%d,%d)",
          a.shape(0), a.shape(1), a.shape(2),
          b.shape(0), b.shape(1), b.shape(2),
          (a.shape(0) == b.shape(0)) && 
          (a.shape(1) == b.shape(1)) &&
          (a.shape(2) == b.shape(2)));

  const int mz = a.shape(0);
  const int my = a.shape(1);
  const int mx = a.shape(2);

  bool all_match = true;
  for (int iz = 0; iz < mz; iz++){
    for (int iy = 0; iy < my; iy++){
      for (int ix = 0; ix < mx; ix++){
        bool match = (a(iz,iy,ix) == b(iz,iy,ix));
        if (!match){
          if (all_match){
            CkPrintf("\nUnequal array element\n");
          }
          CkPrintf("One array has a value of %e at (%d,%d,%d), the other has "
                   "a value of %e",
                   static_cast<double>(a(iz,iy,ix)), iz,iy,ix,
                   static_cast<double>(b(iz,iy,ix)));
          all_match = false;
        }
      }
    }
  }
  Unit::instance()->assertion(all_match, file, line, true);
}

//----------------------------------------------------------------------

template<typename T>
void assert_allequal3D_(const CelloView<T,3>& a, const T& b,
                        const char* file, int line) {

  const int mz = a.shape(0);
  const int my = a.shape(1);
  const int mx = a.shape(2);

  bool all_equal = true;
  for (int iz = 0; iz < mz; iz++){
    for (int iy = 0; iy < my; iy++){
      for (int ix = 0; ix < mx; ix++){
        all_equal &= (a(iz,iy,ix) == b);
      }
    }
  }
  ASSERT("assert_allequal3D", "The elements are not all equal", all_equal);
}

//----------------------------------------------------------------------

template<typename T>
void assert_allequal3D_(const CelloView<T,3>& a, ValRange<T> val_range,
                        const char* file, int line) {

  const int mz = a.shape(0);
  const int my = a.shape(1);
  const int mx = a.shape(2);

  int count = 0;

  bool all_equal = true;
  for (int iz = 0; iz < mz; iz++){
    for (int iy = 0; iy < my; iy++){
      for (int ix = 0; ix < mx; ix++){
        T expected = val_range.start + static_cast<T>(count)*val_range.step;
        count++;

        all_equal &= (a(iz,iy,ix) == expected);
      }
    }
  }
  ASSERT("assert_allequal3D", "The elements are not all equal", all_equal);
}

//----------------------------------------------------------------------

#define assert_allequal3D(VALS, REF)                     \
  assert_allequal3D_(VALS, REF, __FILE__, __LINE__);

//----------------------------------------------------------------------

template<typename T>
void assert_shape3D(const CelloView<T,3>& a,
                    const std::array<int,3>& expected){
  ASSERT("assert_shape3D", "The array doesn't have the expected shape",
         (a.shape(0) == expected[0]) &&
         (a.shape(1) == expected[1]) &&
         (a.shape(2) == expected[2]));
}

#endif /* TEST_VIEW_TEST_ROUNTINES_HPP */
