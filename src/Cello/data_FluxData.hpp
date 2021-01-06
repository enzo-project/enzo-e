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
      if (fd.get_block_fluxes_(i) != nullptr)
        block_fluxes_[i] = new FaceFluxes(*fd.get_block_fluxes_(i));
      if (fd.get_neighbor_fluxes_(i) != nullptr)
        neighbor_fluxes_[i] = new FaceFluxes(*fd.get_neighbor_fluxes_(i));
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
  
  /// Allocate all flux arrays for each field in the list of field
  /// indices.  Optional arrays to indicate the centering of fields
  /// may also be provided
  void allocate
  (int nx, int ny, int nz,
   std::vector<int> field_list,
   std::vector<int> * cx_list=nullptr,
   std::vector<int> * cy_list=nullptr,
   std::vector<int> * cz_list=nullptr);

  /// Deallocate all face fluxes for all faces and all fields
  void deallocate();

  /// Return the number of field indices
  inline unsigned num_fields () const
  { return field_list_.size(); }

  /// Return the i'th field index
  inline int index_field (unsigned i_f) const
  {
    return (0 <= i_f && i_f < field_list_.size()) ?
      field_list_[i_f] : -1;
  }

  /// Return the face fluxes object associated with the given facet
  /// and field.  Note 0 <= i_f < num_fields() is an index into the
  /// field_list vector, not the field index itself.
  inline const FaceFluxes * block_fluxes (int axis, int face, unsigned i_f) const
  {  return get_block_fluxes_ (index_(axis,face,i_f)); }

  inline FaceFluxes * block_fluxes (int axis, int face, unsigned i_f) 
  {  return get_block_fluxes_ (index_(axis,face,i_f)); }

  /// Return the neighboring block's face fluxes associated with the
  /// given facet and field.  Note 0 <= i_f < num_fields() is an index
  /// into the field_list vector, not the field index itself.
  inline const FaceFluxes * neighbor_fluxes (int axis, int face, unsigned i_f) const
  { return get_neighbor_fluxes_ (index_(axis,face,i_f)); }

  inline FaceFluxes * neighbor_fluxes (int axis, int face, unsigned i_f) 
  { return get_neighbor_fluxes_ (index_(axis,face,i_f)); }

  /// Set the block's face fluxes associated with the given facet and
  /// field. Note 0 <= i_f < num_fields() is an index into the
  /// field_list vector, not the field index itself.
  inline void set_block_fluxes
  (FaceFluxes * fluxes, int axis, int face, unsigned i_f)
  { set_block_fluxes(fluxes,index_(axis,face,i_f)); }

  inline void set_block_fluxes (FaceFluxes * fluxes, unsigned i)
  {  block_fluxes_[i] = fluxes; }

  /// Set the neighboring block's face fluxes associated with the
  /// given facet and field. Note 0 <= i_f < num_fields() is an index
  /// into the field_list vector, not the field index itself.
  inline void set_neighbor_fluxes
  (FaceFluxes * fluxes, int axis, int face, unsigned i_f)
  { set_neighbor_fluxes(fluxes,index_(axis,face,i_f)); }

  inline void set_neighbor_fluxes (FaceFluxes * fluxes, unsigned i)
  { neighbor_fluxes_[i] = fluxes; }

  /// Accumulate (sum) the neighboring block's face fluxes associated
  /// with the given facet and field. Note 0 <= i_f < num_fields() is
  /// an index into the field_list vector, not the field index itself.
  inline void sum_neighbor_fluxes
  (FaceFluxes * fluxes, int axis, int face, unsigned i_f)
  { sum_neighbor_fluxes(fluxes,index_(axis,face,i_f)); }

  void sum_neighbor_fluxes (FaceFluxes * fluxes, unsigned index)
  {
    int mx,my,mz;
    neighbor_fluxes_[index]->get_size(&mx,&my,&mz);
    auto & flux_array = fluxes->flux_array();
    auto & neighbor_flux_array = neighbor_fluxes_[index]->flux_array();
    for (int i=0; i<mx*my*mz; i++) {
      neighbor_flux_array[i] += flux_array[i];
    }
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

protected: // functions

  /// Return the block's ith flux data object
  inline const FaceFluxes * get_block_fluxes_ (unsigned i) const
  { return (0 <= i && i < block_fluxes_.size()) ?
      block_fluxes_.at(i) : nullptr; }
  inline FaceFluxes * get_block_fluxes_ (unsigned i) 
  { return (0 <= i && i < block_fluxes_.size()) ?
      block_fluxes_.at(i) : nullptr; }

  /// Return the neighbor block's ith flux data object
  inline const FaceFluxes * get_neighbor_fluxes_ (unsigned i) const
  { return (0 <= i && i < neighbor_fluxes_.size()) ?
      neighbor_fluxes_.at(i) : nullptr;  }
  inline FaceFluxes * get_neighbor_fluxes_ (unsigned i) 
  { return (0 <= i && i < neighbor_fluxes_.size()) ?
      neighbor_fluxes_.at(i) : nullptr;  }

  /// Index for fluxes in block_fluxes_ and neighbor_fluxes_ vectors
  inline int index_ (int axis, int face, unsigned i_f) const
  { return axis + 3*(face + 2*i_f); }

protected: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Face fluxes for this block on each face
  std::vector<FaceFluxes *> block_fluxes_;
  
  /// Face fluxes for neighboring blocks on each face
  std::vector<FaceFluxes *> neighbor_fluxes_;

  /// List of field indices for fluxes
  std::vector<int> field_list_;

};

#endif /* DATA_FLUX_DATA_HPP */

