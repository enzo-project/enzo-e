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
  }

  virtual ~FluxData()
  {
    int n=block_fluxes_.size();
    for (int i=0; i<n; i++) {
      delete block_fluxes_[i];
      delete neighbor_fluxes_[i];
      block_fluxes_[i] = nullptr;
      neighbor_fluxes_[i] = nullptr;
    }
  }

  FluxData( const FluxData & fd )
  {
    int n = fd.block_fluxes_.size();
    block_fluxes_.resize(n);
    neighbor_fluxes_.resize(n);
    // deep-copy fd
    for (int i=0; i<n; i++) {
      if (fd.block_fluxes(i) != nullptr)
        block_fluxes_[i] = new FaceFluxes(*fd.block_fluxes(i));
      if (fd.neighbor_fluxes(i) != nullptr)
        neighbor_fluxes_[i] = new FaceFluxes(*fd.neighbor_fluxes(i));
    }
    field_list_ = fd.field_list_;
    field_index_map_ = fd.field_index_map_;
  }
    
  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change

    const bool pk = p.isPacking();

    int n;
    if (pk) {
      n=block_fluxes_.size();
      for (int i=0; i<n; i++) {
        ASSERT2("FluxData::pup()",
                "block_fluxes_ should be empty but [%d] = %p",
                i,block_fluxes_[i],
                (block_fluxes_[i] == nullptr));
      }
      n=neighbor_fluxes_.size();
      // should be empty
      for (int i=0; i<n; i++) {
        ASSERT2("FluxData::pup()",
                "neighbor_fluxes_ should be empty but [%d] = %p",
                i,neighbor_fluxes_[i],
                (neighbor_fluxes_[i] == nullptr));
      }
    }
    p | field_list_;
    p | field_index_map_;
  }
  
  /// Return the block's ith flux data object
  FaceFluxes * block_fluxes (int i)
  {
    return block_fluxes_[i]; }
  
  const FaceFluxes * block_fluxes (int i) const
  { return block_fluxes_[i]; }

  /// Return the block's fluxes for the given face
  FaceFluxes * block_fluxes (int axis, int face, int index_field)
  {
    auto it = field_index_map_.find(index_field);
    if (it != field_index_map_.end()) {
      const int i_f = field_index_map_[index_field];
      const int i = index_ (axis,face,i_f);
      return block_fluxes_[i];
    } else {
      return nullptr;
    }
  }

  /// Allocate all facet fluxes for all fields
  void allocate
  (int nx, int ny, int nz,
   int level, double dt,
   std::vector<int> field_list,
   std::vector<int> * cx_list=nullptr,
   std::vector<int> * cy_list=nullptr,
   std::vector<int> * cz_list=nullptr);
                 
  /// Return the neighbor block's ith flux data object
  FaceFluxes * neighbor_fluxes (int i)
  { return neighbor_fluxes_[i]; }
  const FaceFluxes * neighbor_fluxes (int i) const
  { return neighbor_fluxes_[i]; }

  /// Return the neighbor block's fluxes for the given face
  FaceFluxes * neighbor_fluxes (int axis, int face, int index_field)
  {
    auto it = field_index_map_.find(index_field);
    if (it != field_index_map_.end()) {
      const int i_f = field_index_map_[index_field];
      const int i = index_ (axis,face,i_f);
      return neighbor_fluxes_[i];
    } else {
      return nullptr;
    }
  }

  /// Set the block's fluxes for the given face
  void set_block_fluxes
  (FaceFluxes * fluxes, int axis, int face, int index_field)
  {
    auto it = field_index_map_.find(index_field);
    if (it != field_index_map_.end()) {
      const int i_f = field_index_map_[index_field];
      const int i = index_ (axis,face,i_f);
      block_fluxes_[i] = fluxes;
    }
  }

  /// Set the block's fluxes for the given face
  void set_neighbor_fluxes
  (FaceFluxes * fluxes, int axis, int face, int index_field)
  {
    auto it = field_index_map_.find(index_field);
    if (it != field_index_map_.end()) {
      const int i_f = field_index_map_[index_field];
      const int i = index_ (axis,face,i_f);
      neighbor_fluxes_[i] = fluxes;
    }
  }

  /// Delete the block's face fluxes object for the given face
  void delete_block_fluxes (int axis, int face, int index_field)
  {
    auto it = field_index_map_.find(index_field);
    if (it != field_index_map_.end()) {
      const int i_f = field_index_map_[index_field];
      const int i = index_ (axis,face,i_f);
      delete block_fluxes_[i];
      block_fluxes_[i] = nullptr;
    }
  }

  /// Delete the block's face fluxes object for the given face
  void delete_neighbor_fluxes (int axis, int face, int index_field)
  {
    auto it = field_index_map_.find(index_field);
    if (it != field_index_map_.end()) {
      const int i_f = field_index_map_[index_field];
      const int i = index_ (axis,face,i_f);
      delete neighbor_fluxes_[i];
      neighbor_fluxes_[i] = nullptr;
    }
  }

  int field_index (int index_field)
  { return field_index_map_[index_field]; }

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

  inline int index_ (int axis, int face, int i_f)
  { return axis + 3*(face + 2*i_f); }
private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Face fluxes for this block on each face
  std::vector<FaceFluxes *> block_fluxes_;
  
  /// Face fluxes for neighboring blocks on each face
  std::vector<FaceFluxes *> neighbor_fluxes_;

  /// List of field indices for fluxes
  std::vector<int> field_list_;

  /// Map of field indices for fluxes (inverse of field_list_)
  std::map<int,int> field_index_map_;

};

#endif /* DATA_FLUX_DATA_HPP */

