#ifndef ENZO_ENZO_PERMUTED_COORDINATES_HPP
#define ENZO_ENZO_PERMUTED_COORDINATES_HPP


// The EnzoPermutedCoordinates represents a Cyclic Permutation of the Cartesian
// Coordinate basis. It has axes i, j, and k.
//
// Allowed permutations are:
//   (k,j,i) = (z,y,x)
//   (k,j,i) = (x,z,y)
//   (k,j,i) = (y,x,z)
//
// The class assumes that each cartesian axis has the following id:
//    0 <--> x, 1 <--> y, 2 <--> z

class EnzoPermutedCoordinates
{

  /// @class    EnzoPermutedCoordinates
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates operations related to the cyclic permutation of the Cartesian Axes

public: // interface

  /// Construct an instance of EnzoPermutedCoordinates by specifying the
  /// axis_id (0, 1, or 2) of the Cartesian axis that lies along axis i
  EnzoPermutedCoordinates(int i_axis)
    : i_axis_(i_axis)
  {
    ASSERT1("EnzoPermutedCoordinates",
	    "Constructor arg must be 0, 1, or 2. Not %d.",
	    i_axis, i_axis>-1 && i_axis<3);
  }

  /// Returns the code of the axis id associated with each direction
  /// 0 <--> x, 1 <--> y, 2 <--> z 
  inline int i_axis() const { return i_axis_; }
  inline int j_axis() const { return (i_axis_+1)%3; }
  inline int k_axis() const { return (i_axis_+2)%3; }

  /// calculates the components of the i, j, and k unit vectors
  inline void i_unit_vector(int &i_z, int &i_y, int &i_x) const {
    unit_vector_comp_(i_axis(), i_z, i_y, i_x);
  }
  inline void j_unit_vector(int &j_z, int &j_y, int &j_x) const {
    unit_vector_comp_(j_axis(), j_z, j_y, j_x);
  }
  inline void k_unit_vector(int &k_z, int &k_y, int &k_x) const {
    unit_vector_comp_(k_axis(), k_z, k_y, k_x);
  }

  /// Returns the subarray of array using the offset start values indicated
  /// using the i,j,k coordinate system represented by this instance
  ///
  /// This method assumes that spatial dimensions of array are orderred as:
  /// (z-axis, y-axis, x-axis)
  ///
  /// Examples
  /// --------
  ///   EnzoPermutedCoordinates(1).left_edge_offset(array, 3, 4, 1) forwards to:
  ///     array.subarray(CSlice(4, array.shape(0)),
  ///                    CSlice(1, array.shape(1)),
  ///                    CSlice(3, array.shape(2)))
  ///
  ///   EnzoPermutedCoordinates(2).left_edge_offset(array, 0, 1, 2) forwards to:
  ///     array.subarray(CSlice(2, array.shape(0)),
  ///                    CSlice(0, array.shape(1)),
  ///                    CSlice(1, array.shape(2)))
  template<typename T>
  inline CelloView<T,3> left_edge_offset(const CelloView<T,3> &array,
                                         int kstart, int jstart, int istart)
    const noexcept
  {
    return get_subarray(array, CSlice(kstart, nullptr), CSlice(jstart, nullptr),
                        CSlice(istart,nullptr));
  }

  /// Returns the subarray of array for slices specified along the k-, j-,
  /// and i- axes using the i,j,k coordinate system represented by this instance
  template<typename T>
  inline CelloView<T,3> get_subarray(const CelloView<T,3> &array,
                                     CSlice k_slice, CSlice j_slice,
                                     CSlice i_slice) const noexcept
  {
    CSlice slices[3];
    slices[2 - i_axis()] = i_slice;
    slices[2 - j_axis()] = j_slice;
    slices[2 - k_axis()] = k_slice;
    return array.subarray(slices[0], slices[1], slices[2]);
  }

private:

  // helper_function
  inline void unit_vector_comp_(int v_axis, int &v_k, int &v_j, int &v_i) const
  {
    v_i = (v_axis == 0) ? 1 : 0;
    v_j = (v_axis == 1) ? 1 : 0;
    v_k = (v_axis == 2) ? 1 : 0;
  }

private: // attributes
  int i_axis_;
};

#endif /*ENZO_ENZO_PERMUTED_COORDINATES_HPP*/
