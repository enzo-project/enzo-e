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

  /// Create an empty FluxData() object

  FluxData()
    : face_flux_map_()
  { }

  /// Insert the given FaceFluxes object into the FluxData. Must be
  /// dynamically allocated, and responsibility for deleting
  /// face_fluxes is transfered to the FluxData object.

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change

    //p | face_flux_map_;

  }

  /// Insert the FaceFluxes object indexed by face
  
  void insert_fluxes(Face face, FaceFluxes * face_fluxes)
  {
    auto it = face_flux_map_.find(face);
    if (it == face_flux_map_.end()) {
      face_flux_map_.insert
        (std::pair<Face,FaceFluxes*>(face,face_fluxes));
    }
  }

  /// Remove the FaceFluxes object indexd by face from the FluxData object,
  /// but do not delete it.

  void remove_fluxes(Face face, FaceFluxes * face_fluxes)
  {
    auto it = face_flux_map_.find(face);
    if (it != face_flux_map_.end()) {
      face_flux_map_.erase(it);
    }
  }

  /// Remove the FaceFluxes object indexed by face from the FluxData
  /// object and delete it.

  void delete_fluxes(Face face, FaceFluxes * face_fluxes)
  {
    auto it = face_flux_map_.find(face);
    if (it != face_flux_map_.end()) {
      face_flux_map_.erase(it);
      delete it->second;
    }
  }

  /// Return the FaceFluxes object for the given Face, Block index,
  /// and Field index. Block index must match one of defining Block
  /// indices for Face.

  FaceFluxes face_fluxes (Face face, Index index_block, int index_field)
  {
    auto it = face_flux_map_.find(face);
    if (it != face_flux_map_.end()) {
      return *it->second;
    }
  }

  /// Return the number of FaceFluxes in the FluxData object.
  size_t num_face_fluxes()
  {
    return face_flux_map_.size();
  }

  /// Return the ith FaceFluxes object.
  FaceFluxes * face_fluxes (int n)
  {
    std::map<Face, FaceFluxes * >::iterator it;
    for (it=face_flux_map_.begin() ;
           (it != face_flux_map_.end()) && --n>=0;
           it++)
      ;
    return it == face_flux_map_.end() ? nullptr : it->second;
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

  /// FaceFluxes indexed by Face objects
  std::map<Face, FaceFluxes * > face_flux_map_;
};

#endif /* DATA_FLUX_DATA_HPP */

