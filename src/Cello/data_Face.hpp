// See LICENSE_CELLO file for license and copyright information

/// @file     data_Face.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2019-10-08
/// @brief    [\ref Data] Declaration of the Face class

#ifndef DATA_FACE_HPP
#define DATA_FACE_HPP

class Face {

  /// @class    Face
  /// @ingroup  Data
  /// @brief    [\ref Data]

public: // interface

  Face ()
  { }
  
  /// Constructor
  Face(int ix, int iy, int iz, int axis, int face) throw()
    :  ix_(ix),iy_(iy),iz_(iz),
       axis_(axis), face_(face)
  {
  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;

    p | ix_;
    p | iy_;
    p | iz_;
    p | axis_;
    p | face_;
  }

  /// Return the tuple (ix,iy,iz), -1 <= ix,iy,iz <= 1, identifying
  /// the block's face, which may be a corner, edge, facet, or the
  /// entire block
  void get_face (int *ix, int *iy, int *iz) const
  {
    if (ix) (*ix) = ix_;
    if (iy) (*iy) = iy_;
    if (iz) (*iz) = iz_;
  }

  /// Return the axis associated with the normal direction: x=0, y=1,
  /// or z=2.*
  inline int axis() const
  { return axis_; }
  
  /// Return whether the normal direction is towards the lower (0)
  /// or upper (1) face direction.*
  inline int face() const
  { return face_; }
  
  //--------------------------------------------------

  /// Return the number of bytes required to serialize the data object
  int data_size () const;

  /// Serialize the object into the provided empty memory buffer.
  /// Returns the next open position in the buffer to simpliiy
  /// serializing multiple objects in one buffer.
  char * save_data (char * buffer) const;

  /// Restore the object from the provided initialized memory buffer data.
  /// Returns the next open position in the buffer to simpliiy
  /// serializing multiple objects in one buffer.
  char * load_data (char * buffer);

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Face (facet, edge, or corner) specified by -1 <= ix_,iy_,iz_ <= 1
  int ix_,iy_,iz_;

  /// Axis 0 <= axis_ < rank associated with the normal vector
  int axis_;

  /// Face 0 (lower) or 1 (upper) associated with the normal vector
  int face_;
};

#endif /* DATA_FACE_HPP */

