// See LICENSE_CELLO file for license and copyright information

/// @file     test_ArrayTestTools.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Fri Jan 20 2023
/// @brief    [\ref Test] Define helper functions to assist with unit testing
///
/// These functions primarily assist with tests involving CelloViews

//----------------------------------------------------------------------

/// fill an existing CelloView with values that start at ``start`` and are
/// incremented by ``step``
template<typename T>
void assign_range_3Darray(const CelloView<T,3>& arr, T start, T step) {
  int count = 0;

  const int mz = arr.shape(0);
  const int my = arr.shape(1);
  const int mx = arr.shape(2);

  for (int iz = 0; iz < mz; iz++){
    for (int iy = 0; iy < my; iy++){
      for (int ix = 0; ix < mx; ix++){
        arr(iz,iy,ix) = start + static_cast<T>(count)*step;
        count++;
      }
    }
  }

}

//----------------------------------------------------------------------

/// allocate a new CelloView and fill it with values that start from ``start``
/// and are incremented by ``step``
template<typename T>
CelloView<T,3> range_3Darray(int mz, int my, int mx, T start, T step) {
  CelloView<T,3> out(mz,my,mx);
  assign_range_3Darray(out, start, step);
  return out;
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

#define assert_allequal3D(VALS, REF)                     \
  assert_allequal3D_(VALS, REF, __FILE__, __LINE__);

