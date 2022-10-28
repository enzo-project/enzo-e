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
   int index_refine,
   int num_refine);

  EnzoMethodInference()
    : Method(),
      level_base_(0),
      level_array_(0),
      level_infer_(0),
      m3_level_(),
      m3_infer_(),
      field_group_(),
      is_sync_child_(-1),
      is_sync_parent_(-1),
      is_mask_(-1),
      index_refine_(-1),
      num_refine_(0)
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
      is_sync_child_(-1),
      is_sync_parent_(-1),
      is_mask_(-1),
      index_refine_(-1),
      num_refine_(0)
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
  void create_arrays (Block * block);

protected: // methods

  /// Apply criteria to determine which if any overlapping inference
  /// arrays need to be created. Return number of inference arrays
  /// to create, and allocate and initialize inference mask
  int apply_criteria_(Block * block);

  /// Allocate the inference array mask for the given block if not already
  /// allocated, and return the mask array length
  int mask_allocate_(Block * block);
  /// Get the mask size for blocks in the given level, and return
  /// the mask array length
  int mask_dims_(int level,int *mx=nullptr, int *my=nullptr, int *mz=nullptr) const;

  /// Create the level arrays according to the Block's level array mask,
  /// and return the number of level arrays created
  void create_level_arrays_ (Block * block);

  /// Process inference array creation criteria results, sending
  ///results where needed so that inference arrays are created
  void forward_create_array_ (Block * block, int count);

  void compute_ (Block * block, enzo_float * Unew ) throw();

  void intersecting_root_blocks_
  (int ib3_lower[3], int ib3_upper[3], const int ia3[3],
   const int n3[3], int level,
   const int na3[3], const int nb3[3]) const;

  /// Return the dimensionality of the level array
  void level_array_dims_(int *mx, int *my, int *mz);
 
  void intersecting_level_arrays_
  (int ia3_lower[3], int ia3_upper[3], const int ib3[3],
   const int n3[3], int level,
   const int na3[3], const int nb3[3]) const;

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

  /// Print the char scalar array mask
  void print_mask_(Block * block) const
  {
    int mx,my,mz;
    int n = mask_dims_(block->level(),&mx,&my,&mz);
    const char * mask = *scalar_mask_(block);
    for (int iy=my-1; iy>=0; iy--) {
      CkPrintf ("MASK %s ",block->name().c_str());
      for (int ix=0; ix<mx; ix++) {
        int c=0;
        for (int iz=0; iz<mz; iz++) {
          int i=ix+mx*(iy+my*iz);
          c+= (mask[i])?1:0;
        }
        CkPrintf ("%1d",c);
      }
      CkPrintf ("\n");
    }
  }
  
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

  /// Block counter index for child synchronization
  int is_sync_child_;

  /// Block counter index for parent synchronization
  int is_sync_parent_;

  /// Block Scalar mask index for creating inference arrays
  int is_mask_;

  /// Index for the first "refinement" criterion
  int index_refine_;
  /// Number of refinement criteria
  int num_refine_;
};

#endif /* ENZO_ENZO_METHOD_INFERENCE_HPP */
