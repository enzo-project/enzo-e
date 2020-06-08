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
  Face(int ix, int iy, int iz,
       int rx, int ry, int rz) throw()
    :  ix_(ix),iy_(iy),iz_(iz),
       rx_(rx),ry_(ry),rz_(rz)
  {
  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;

    p | ix_;
    p | iy_;
    p | iz_;
    p | rx_;
    p | ry_;
    p | rz_;
  }

  inline bool operator==(const Face& f)
  {
    return (ix_ == f.ix_)
      &&   (iy_ == f.iy_)
      &&   (iz_ == f.iz_)
      &&   (rx_ == f.rx_)
      &&   (ry_ == f.ry_)
      &&   (rz_ == f.rz_);
  }

  void get_face (int *ix, int *iy, int *iz) const
  {
    if (ix) (*ix) = ix_;
    if (iy) (*iy) = iy_;
    if (iz) (*iz) = iz_;
  }

  void get_normal (int *rx, int *ry, int *rz) const
  {
    if (rx) (*rx) = rx_;
    if (ry) (*ry) = ry_;
    if (rz) (*rz) = rz_;
  }

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
  /// Normal direction
  int rx_,ry_,rz_;
};

#endif /* DATA_FACE_HPP */

