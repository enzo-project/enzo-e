// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodInference.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Thu Apr  1 16:14:38 PDT 2010
/// @brief    [\ref Enzo] Declaration of EnzoMethodInference
///           forward Euler solver for the Inference equation

#ifndef ENZO_ENZO_METHOD_INFERENCE_HPP
#define ENZO_ENZO_METHOD_INFERENCE_HPP

class EnzoMethodInference : public Method {

  /// @class    EnzoMethodInference
  /// @ingroup  Enzo
  ///
  /// @brief [\ref Enzo] Demonstration method to solve Inference equation
  /// using forward Euler method

public: // interface

  /// Create a new EnzoMethodInference object
  EnzoMethodInference
  (int level_base,
   int level_array,
   int level_infer,
   const int array_dims[3],
   const int array_size[3],
   std::string field_group,
   float overdensity_threshold);

  EnzoMethodInference()
    : Method(),
      level_base_(0),
      level_array_(0),
      level_infer_(0),
      m3_level_(),
      m3_infer_(),
      field_group_(),
      num_fields_(0),
      is_sync_child_(-1),
      is_sync_parent_(-1),
      is_mask_(-1),
      is_count_(-1),
      overdensity_threshold_(0)
  { }

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodInference);

  /// Charm++ PUP::able migration constructor
  EnzoMethodInference (CkMigrateMessage *m)
    : Method (m),
      level_base_(0),
      level_array_(0),
      level_infer_(0),
      m3_level_(),
      m3_infer_(),
      field_group_(),
      num_fields_(0),
      is_sync_child_(-1),
      is_sync_parent_(-1),
      is_mask_(-1),
      is_count_(-1),
      overdensity_threshold_(0)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) ;

  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "inference"; }

public: // methods

  void merge_masks (Block * block, int n, char *mask, int ic3[3]);

  void count_arrays (Block * block, int count);
  // void create_arrays (Block * block);

  /// Return the field group defining fields in the inference arrays
  std::string field_group () const
  { return field_group_; }

  /// Request from inference array index ia3[3] for data from given
  /// block
  void request_data (Block * block, int ia3[3]);

  void update (Block *, int n, char * buffer, int ia3[3]);

protected: // methods

  /// Apply criteria to determine which if any overlapping inference
  /// arrays need to be created, and create an associated inference
  /// array mask
  void apply_criteria_(Block * block);

  /// Allocate the inference array mask for the given block if not already
  /// allocated
  void mask_allocate_(Block * block, int nx, int ny, int nz);

  /// Get the mask size for blocks in the given level, and return
  /// the mask array length
  std::tuple<int,int,int> mask_dims_(int level) const;

  /// Compute local overdensity, tagging mask array accordingly
  void compute_overdensity_
  (Block * block, char * mask, int nx, int ny, int nz);

  /// Compute offsets ox,oy,oz and sizes nx,ny,nz into field
  /// data, taking into account one ghost zone layer, field
  /// ghost zone depth, restriction count, block size relative
  /// to inference array, etc.
  std::tuple<int,int,int, int,int,int>
  get_block_portion_
  (Index index, int index_field, int ia3[3]);

/// Create the level arrays according to the Block's level array mask,
  /// and return the number of level arrays created
  void create_level_arrays_ (Block * block);

  /// Process inference array creation criteria results, sending
  ///results where needed so that inference arrays are created
  void forward_create_array_ (Block * block, int count);

  /// Merge level array masks of children into this block
  void merge_masks_1to1_ (char * mask, char * mask_in, int level);
  void merge_masks_2to1_ (char * mask, char * mask_in, int level, int ic3[3]);

  /// Count inference arrays to create given block array mask
  int count_arrays_ (Block * block) const;

  void compute_ (Block * block, enzo_float * Unew ) throw();

  bool block_intersects_array_(Index index, int ia3[3]);

  /// Return the dimensionality of the level array
  void level_array_dims_(int *mx, int *my, int *mz);
 
  /// Return the Block's synchronization counter for child blocks (plus self)
  Sync & sync_child_(Block * block)
  { return *block->data()->scalar_sync().value(is_sync_child_); }

  /// Return the Block's synchronization counter for parent block (plus self)
  Sync & sync_parent_(Block * block)
  { return *block->data()->scalar_sync().value(is_sync_parent_); }

  /// Return a pointer to the char scalar array mask
  char ** scalar_mask_(Block * block)
  { return (char **)block->data()->scalar_void().value(is_mask_); }
  const char ** scalar_mask_(Block * block) const
  { return (const char **)block->data()->scalar_void().value(is_mask_); }

  /// Return count of overlapping inference arrays in a block
  int & scalar_count_(Block * block)
  { return *block->data()->scalar_int().value(is_count_); }

  void coarsen_
  (enzo_float * ac,
   int mcx, int mcy, int mcz, int ncx, int ncy, int ncz, int ecx, int ecy, int ecz,
   const enzo_float * af,
   int mfx, int mfy, int mfz, int nfx, int nfy, int nfz, int efx, int efy, int efz);

protected: // attributes

  /// Base level of blocks interacting with the level array
  int level_base_;

  /// Level at which the level array lives (one element per block at
  /// this level)
  int level_array_;

  /// Level of the resolution for the inference arrays.
  int level_infer_;

  /// Dimensions of the level array, each element of which contains
  /// inferenece arrays for each required field
  int m3_level_[3];

  /// Size of the inference arrays associated with each level_array element
  int m3_infer_[3];

  /// Field group defining which fields to include in the inference
  /// array
  std::string field_group_;

  /// Number of fields in the field group
  int num_fields_;

  /// Block counter index for child synchronization
  int is_sync_child_;

  /// Block counter index for parent synchronization
  int is_sync_parent_;

  /// Block Scalar mask index for creating inference arrays
  int is_mask_;

  /// Block Scalar count of overlapping inference arrays
  int is_count_;

  /// Local overdensity threshold for creating inference array
  enzo_float overdensity_threshold_;

};

#endif /* ENZO_ENZO_METHOD_INFERENCE_HPP */
