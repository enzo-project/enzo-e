// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodInferenceArray.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Thu Apr  1 16:14:38 PDT 2010
/// @brief    [\ref Enzo] Declaration of EnzoMethodInferenceArray
///           forward Euler solver for the InferenceArray equation

#ifndef ENZO_ENZO_METHOD_INFERENCE_ARRAY_HPP
#define ENZO_ENZO_METHOD_INFERENCE_ARRAY_HPP

class EnzoMethodInferenceArray : public Method {

  /// @class    EnzoMethodInferenceArray
  /// @ingroup  Enzo
  ///
  /// @brief [\ref Enzo] Demonstration method to solve InferenceArray equation
  /// using forward Euler method

public: // interface

  /// Create a new EnzoMethodInferenceArray object
  EnzoMethodInferenceArray
  (int level,
   const int root_size[3],
   const int array_dims[3],
   const int array_size[3],
   const int array_ghosts[3],
   std::string field_group);

  EnzoMethodInferenceArray()
    : Method(),
      level_(0),
      array_dims_(),
      array_size_(),
      array_ghosts_(),
      field_group_(),
      is_sync_child_(-1),
      is_sync_parent_(-1)
  { }

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodInferenceArray);

  /// Charm++ PUP::able migration constructor
  EnzoMethodInferenceArray (CkMigrateMessage *m)
    : Method (m),
      level_(0),
      array_dims_(),
      array_size_(),
      array_ghosts_(),
      field_group_(),
      is_sync_child_(-1),
      is_sync_parent_(-1)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) ;

  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "inference_array"; }

protected: // methods

  void compute_ (Block * block, enzo_float * Unew ) throw();

  void intersecting_root_blocks_
  (int ib3_lower[3], int ib3_upper[3], const int ia3[3],
   const int n3[3], int level,
   const int na3[3], const int ga3[3], const int nb3[3]) const;

  void intersecting_level_arrays_
  (int ia3_lower[3], int ia3_upper[3], const int ib3[3],
   const int n3[3], int level,
   const int na3[3], const int ga3[3], const int nb3[3]) const;

  /// Return the Block's synchronization counter for child blocks (plus self)
  Sync & s_sync_child_(EnzoBlock * block)
  { return *block->data()->scalar_sync().value(is_sync_child_); }

  /// Return the Block's synchronization counter for parent block (plus self)
  Sync & s_sync_parent_(EnzoBlock * block)
  { return *block->data()->scalar_sync().value(is_sync_parent_); }

protected: // attributes

  /// Refinement level corresponding to the resolution of the inference array
  int level_;

  /// Dimensions of the level array, each element contains inferenece arrays
  int array_dims_[3];

  /// Size of the inference arrays associated with each level_array element
  int array_size_[3];

  /// Number of overlapping level_ ghost cells in the inference arrays
  int array_ghosts_[3];

  /// Field group defining which fields to include in the inference
  /// array
  std::string field_group_;

  /// Block counter for child synchronization
  int is_sync_child_;

  /// Block counter for parent synchronization
  int is_sync_parent_;
};

#endif /* ENZO_ENZO_METHOD_INFERENCE_ARRAY_HPP */
