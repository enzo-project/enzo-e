// See LICENSE_CELLO file for license and copyright information

/// @file     data_FluxData.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2019-10-16
/// @brief    [\ref Data] Declaration of the FluxData class

#ifndef DATA_FLUX_DATA_HPP
#define DATA_FLUX_DATA_HPP

class FluxData {

  /// @class    FluxData
  /// @ingroup  Data
  /// @brief    [\ref Data] 

public: // interface

  //  friend FaceFluxes;
  /// Create an empty FluxData() object

  FluxData()
  {
    
    face_fluxes_block_.resize(27);
    face_fluxes_neighbor_.resize(27);
    std::fill_n(&face_fluxes_block_[0],27,nullptr);
    std::fill_n(&face_fluxes_neighbor_[0],27,nullptr);
  }

  virtual ~FluxData()
  {
    for (int i=0; i<27; i++) {
      delete face_fluxes_block_[i];
      delete face_fluxes_neighbor_[i];
      face_fluxes_block_[i] = nullptr;
      face_fluxes_neighbor_[i] = nullptr;
    }
  }

  FluxData( const FluxData & fd )
  {
    face_fluxes_block_.resize(27);
    face_fluxes_neighbor_.resize(27);
    // deep-copy fd
    for (int i=0; i<27; i++) {
      if (fd.fluxes_block(i) != nullptr)
        face_fluxes_block_[i] = new FaceFluxes(*fd.fluxes_block(i));
      if (fd.fluxes_neighbor(i) != nullptr)
        face_fluxes_neighbor_[i] = new FaceFluxes(*fd.fluxes_neighbor(i));
    }
  }
    
  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change

    const bool pk = p.isPacking();

    int n;
    if (pk) {
      n=face_fluxes_block_.size();
      for (int i=0; i<n; i++) {
        ASSERT2("FluxData::pup()",
                "face_fluxes_block_ should be empty but [%d] = %p",
                i,face_fluxes_block_[i],
                (face_fluxes_block_[i] == nullptr));
      }
      n=face_fluxes_neighbor_.size();
      // should be empty
      for (int i=0; i<n; i++) {
        ASSERT2("FluxData::pup()",
                "face_fluxes_neighbor_ should be empty but [%d] = %p",
                i,face_fluxes_neighbor_[i],
                (face_fluxes_neighbor_[i] == nullptr));
      }
    }
  }
  /// Return the block's ith flux data object
  FaceFluxes * fluxes_block (int i)
  { return face_fluxes_block_[i]; }
  const FaceFluxes * fluxes_block (int i) const
  { return face_fluxes_block_[i]; }

  /// Return the block's fluxes for the given face
  FaceFluxes * fluxes_block (int ix, int iy, int iz)
  {
    const int i = ++ix+3*(++iy+3*++iz);
    return face_fluxes_block_[i];
  }

  /// Return the neighbor block's ith flux data object
  FaceFluxes * fluxes_neighbor (int i)
  { return face_fluxes_neighbor_[i]; }
  const FaceFluxes * fluxes_neighbor (int i) const
  { return face_fluxes_neighbor_[i]; }

  /// Return the neighbor block's fluxes for the given face
  FaceFluxes * fluxes_neighbor (int ix, int iy, int iz)
  {
    const int i = ++ix+3*(++iy+3*++iz);
    return face_fluxes_neighbor_[i]; 
  }

  /// Set the block's fluxes for the given face
  void set_fluxes_block (FaceFluxes * fluxes, int ix, int iy, int iz)
  {
    const int i = ++ix+3*(++iy+3*++iz);
    delete face_fluxes_block_[i];
    face_fluxes_block_[i] = fluxes;
  }

  /// Set the block's fluxes for the given face
  void set_fluxes_neighbor (FaceFluxes * fluxes, int ix, int iy, int iz)
  {
    const int i = ++ix+3*(++iy+3*++iz);
    delete face_fluxes_neighbor_[i];
    face_fluxes_neighbor_[i] = fluxes; 
  }

  /// Delete the block's face fluxes object for the given face
  void delete_fluxes_block (int ix, int iy, int iz)
  {
    const int i = ++ix+3*(++iy+3*++iz);
    delete face_fluxes_block_[i];
    face_fluxes_block_[i] = nullptr;
  }

  /// Delete the block's face fluxes object for the given face
  void delete_fluxes_neighbor (int ix, int iy, int iz)
  {
    const int i = ++ix+3*(++iy+3*++iz);
    delete face_fluxes_neighbor_[i];
    face_fluxes_neighbor_[i] = nullptr;
  }

  //--------------------------------------------------

  /// Return the number of bytes required to serialize the data object
  int data_size () const;

  /// Serialize the object into the provided empty memory buffer.
  /// Returns the next open position in the buffer to simplify
  /// serializing multiple objects in one buffer.
  char * save_data (char * buffer) const;

  /// Restore the object from the provided initialized memory buffer data.
  /// Returns the next open position in the buffer to simplify
  /// serializing multiple objects in one buffer.
  char * load_data (char * buffer);

  //--------------------------------------------------

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Face fluxes for this block on each face
  std::vector<FaceFluxes *> face_fluxes_block_;
  
  /// Face fluxes for neighboring blocks on each face
  std::vector<FaceFluxes *> face_fluxes_neighbor_;
};

#endif /* DATA_FLUX_DATA_HPP */

