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
      index_refine_(-1),
      num_refine_(0)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) ;

  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "inference"; }

protected: // methods

  /// Apply criteria to determine which if any overlapping inference
  /// arrays need to be created. Return number of inference arrays
  /// to create
  int apply_criteria_(Block * block, std::vector<bool> & mask);

  /// Process inference array creation criteria results, sending
  ///results where needed so that inference arrays are created
  void forward_create_array_
  (Block * block, int count, const std::vector<bool> & mask);

  void compute_ (Block * block, enzo_float * Unew ) throw();

  void intersecting_root_blocks_
  (int ib3_lower[3], int ib3_upper[3], const int ia3[3],
   const int n3[3], int level,
   const int na3[3], const int nb3[3]) const;

  void intersecting_level_arrays_
  (int ia3_lower[3], int ia3_upper[3], const int ib3[3],
   const int n3[3], int level,
   const int na3[3], const int nb3[3]) const;

  /// Return the Block's synchronization counter for child blocks (plus self)
  Sync & s_sync_child_(EnzoBlock * block)
  { return *block->data()->scalar_sync().value(is_sync_child_); }

  /// Return the Block's synchronization counter for parent block (plus self)
  Sync & s_sync_parent_(EnzoBlock * block)
  { return *block->data()->scalar_sync().value(is_sync_parent_); }

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

  /// Block counter for child synchronization
  int is_sync_child_;

  /// Block counter for parent synchronization
  int is_sync_parent_;

  /// Index for the first "refinement" criterion
  int index_refine_;
  /// Number of refinement criteria
  int num_refine_;
};

#endif /* ENZO_ENZO_METHOD_INFERENCE_HPP */
