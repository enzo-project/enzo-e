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
    : block_fluxes_(),
      neighbor_fluxes_(),
      field_list_()
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
  }
  
  /// Allocate all facet fluxes for all fields
  void allocate
  (int nx, int ny, int nz,
   int level, double dt,
   std::vector<int> field_list,
   std::vector<int> * cx_list=nullptr,
   std::vector<int> * cy_list=nullptr,
   std::vector<int> * cz_list=nullptr);

  /// Deallocate all face fluxes for all fields
  void deallocate();

  /// Return the block's ith flux data object
  const FaceFluxes * block_fluxes (int axis, int face, int i_f) const
  {  return block_fluxes(index_ (axis,face,i_f)); }
  FaceFluxes * block_fluxes (int axis, int face, int i_f) 
  {  return block_fluxes(index_ (axis,face,i_f)); }
  /// Return the block's fluxes for the given face
  const FaceFluxes * block_fluxes (int i) const
  { return (0 <= i && i < block_fluxes_.size()) ?
      block_fluxes_.at(i) : nullptr; }
  FaceFluxes * block_fluxes (int i) 
  { return (0 <= i && i < block_fluxes_.size()) ?
      block_fluxes_.at(i) : nullptr; }
  /// Return the neighbor block's fluxes for the given face
  const FaceFluxes * neighbor_fluxes (int axis, int face, int i_f) const
  { return neighbor_fluxes(index_(axis,face,i_f)); }
  FaceFluxes * neighbor_fluxes (int axis, int face, int i_f) 
  { return neighbor_fluxes(index_(axis,face,i_f)); }
  /// Return the neighbor block's ith flux data object
  const FaceFluxes * neighbor_fluxes (int i) const
  { return (0 <= i && i < neighbor_fluxes_.size()) ?
      neighbor_fluxes_.at(i) : nullptr;  }
  FaceFluxes * neighbor_fluxes (int i) 
  { return (0 <= i && i < neighbor_fluxes_.size()) ?
      neighbor_fluxes_.at(i) : nullptr;  }

  /// Set the block's fluxes for the given face
  void set_block_fluxes (FaceFluxes * fluxes, int axis, int face, int i_f)
  { set_block_fluxes(fluxes,index_(axis,face,i_f)); }
  void set_block_fluxes (FaceFluxes * fluxes, int i)
  {
    ASSERT2 ("FluxData::set_block_fluxes",
             "Trying to assign to block_fluxes_[%d] vector of size %ud",
             i,block_fluxes_.size(),
             (0 <= i && i < block_fluxes_.size()));
    block_fluxes_[i] = fluxes;
  }

  /// Add the block's fluxes for the given face
  void add_block_fluxes (FaceFluxes * fluxes, int axis, int face, int i_f)
  { add_block_fluxes(fluxes,index_(axis,face,i_f)); }
  void add_block_fluxes (FaceFluxes * fluxes, int index)
  {
    ASSERT2 ("FluxData::add_block_fluxes",
             "Trying to assign to block_fluxes_[%d] vector of size %ud",
             index,block_fluxes_.size(),
             (0 <= index && index < block_fluxes_.size()));
    ASSERT1 ("FluxData::add_block_fluxes",
             "Trying to add to unallocated block_fluxes_[%d]",
             index,
             (block_fluxes_[index] != nullptr));
    int mx,my,mz;
    block_fluxes_[index]->get_dimensions(&mx,&my,&mz);
    auto & flux_array       = fluxes->flux_array();
    auto & block_flux_array = block_fluxes_[index]->flux_array();
    for (int i=0; i<mx*my*mz; i++) {
      block_flux_array[i] += flux_array[i];
    }
  }

  /// Set the neighbor block's fluxes for the given face
  void set_neighbor_fluxes (FaceFluxes * fluxes, int axis, int face, int i_f)
  { set_neighbor_fluxes(fluxes,index_(axis,face,i_f)); }
  void set_neighbor_fluxes (FaceFluxes * fluxes, int i)
  {
    ASSERT2 ("FluxData::set_neighbor_fluxes",
             "Trying to assign to neighbor_fluxes_[%d] vector of size %ud",
             i,neighbor_fluxes_.size(),
             (0 <= i && i < neighbor_fluxes_.size()));
    neighbor_fluxes_[i] = fluxes;
  }

  /// Add the neighborblock's fluxes for the given face
  void add_neighbor_fluxes (FaceFluxes * fluxes, int axis, int face, int i_f)
  { add_neighbor_fluxes(fluxes,index_(axis,face,i_f)); }
  void add_neighbor_fluxes (FaceFluxes * fluxes, int index)
  {
    ASSERT2 ("FluxData::add_neighbor_fluxes",
             "Trying to assign to neighbor_fluxes_[%d] vector of size %ud",
             index,neighbor_fluxes_.size(),
             (0 <= index && index < neighbor_fluxes_.size()));
    ASSERT1 ("FluxData::add_neighbor_fluxes",
             "Trying to add to unallocated neighbor_fluxes_[%d]",
             index,
             (neighbor_fluxes_[index] != nullptr));
    int mx,my,mz;
    neighbor_fluxes_[index]->get_dimensions(&mx,&my,&mz);
    auto & flux_array = fluxes->flux_array();
    auto & neighbor_flux_array = neighbor_fluxes_[index]->flux_array();
    for (int i=0; i<mx*my*mz; i++) {
      neighbor_flux_array[i] += flux_array[i];
    }
  }

  /// Delete the block's face fluxes object for the given face
  void delete_block_fluxes (int axis, int face, int i_f)
  {
    const int i = index_ (axis,face,i_f);
    delete block_fluxes(i);
    block_fluxes_[i] = nullptr;
  }

  /// Delete the block's face fluxes object for the given face
  void delete_neighbor_fluxes (int axis, int face, int i_f)
  {
    const int i = index_ (axis,face,i_f);
    delete neighbor_fluxes(i);
    neighbor_fluxes_[i] = nullptr;
  }

  int index_field (int i_f) const
  {
    return (0 <= i_f && i_f < field_list_.size()) ?
      field_list_[i_f] : -1; }

  int num_fields () const
  { return field_list_.size(); }

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

  /// Index for fluxes in block_fluxes_ and neighbor_fluxes_ vectors
  inline int index_ (int axis, int face, int i_f) const
  { return axis + 3*(face + 2*i_f); }
private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Face fluxes for this block on each face
  std::vector<FaceFluxes *> block_fluxes_;
  
  /// Face fluxes for neighboring blocks on each face
  std::vector<FaceFluxes *> neighbor_fluxes_;

  /// List of field indices for fluxes
  std::vector<int> field_list_;

};

#endif /* DATA_FLUX_DATA_HPP */

