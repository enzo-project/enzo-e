#ifndef ENZO_ENZO_PERMUTED_COORDINATES_HPP
#define ENZO_ENZO_PERMUTED_COORDINATES_HPP


// The EnzoPermutedCoordinates represents a Cyclic Permutation of the Cartesian
// Coordinate basis. It has axes i, j, and k.
//
// Allowed permutations are:
//   (i,j,k) = (x,y,z)
//   (i,j,k) = (y,z,x)
//   (i,j,k) = (z,x,y)
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
	    i_axis>-1 && i_axis<3, i_axis);
  }

  // Returns the code of the axis id associated with each direction
  // 0 <--> x, 1 <--> y, 2 <--> z 
  int i_axis() const { return i_axis_; }
  int j_axis() const { return (i_axis_+1)%3; }
  int k_axis() const { return (i_axis_+2)%3; }

  // calculates the components of the i, j, and k unit vectors
  void i_unit_vector(int &i_x, int &i_y, int &i_z) const {
    unit_vector_comp_(i_axis(), i_x, i_y, i_z);
  }
  void j_unit_vector(int &j_x, int &j_y, int &j_z) const {
    unit_vector_comp_(j_axis(), j_x, j_y, j_z);
  }
  void k_unit_vector(int &k_x, int &k_y, int &k_z) const {
    unit_vector_comp_(k_axis(), k_x, k_y, k_z);
  }

  // Returns the subarray of array using the offset start values indicated
  // using the i,j,k coordinate system represented by this instance
  //
  // This method assumes that spatial dimensions of array are orderred as:
  // (z-axis, y-axis, x-axis)
  //
  // Examples
  // --------
  //   EnzoPermutedCoordinates(1).left_edge_offset(array, 3, 4, 1) forwards to:
  //     array.subarray(ESlice(4, array.shape(0)),
  //                    ESlice(1, array.shape(1)),
  //                    ESlice(3, array.shape(2)))
  //
  //   EnzoPermutedCoordinates(2).left_edge_offset(array, 0, 1, 2) forwards to:
  //     array.subarray(ESlice(2, array.shape(0)),
  //                    ESlice(0, array.shape(1)),
  //                    ESlice(1, array.shape(2)))
  EFlt3DArray left_edge_offset(EFlt3DArray &array, int kstart, int jstart,
			       int istart) const
  {
    int xstart, ystart, zstart;
    if (i_axis() == 0){
      xstart = istart;
      ystart = jstart;
      zstart = kstart;
    } else if (i_axis() == 1){
      ystart = istart;
      zstart = jstart;
      xstart = kstart;
    } else {
      zstart = istart;
      xstart = jstart;
      ystart = kstart;
    }
    return array.subarray(ESlice(zstart, array.shape(0)),
			  ESlice(ystart, array.shape(1)),
			  ESlice(xstart, array.shape(2)));
  }

private:

  // helper_function
  void unit_vector_comp_(int v_axis, int &v_i, int &v_j, int &v_k) const
  {
    v_i = (v_axis == 0) ? 1 : 0;
    v_j = (v_axis == 1) ? 1 : 0;
    v_k = (v_axis == 2) ? 1 : 0;
  }

private: // attributes
  int i_axis_;
};

#endif /*ENZO_ENZO_PERMUTED_COORDINATES_HPP*/
