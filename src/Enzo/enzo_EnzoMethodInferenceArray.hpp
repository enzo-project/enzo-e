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
  (int level, int array_size,int array_ghosts,std::string field_group);

  EnzoMethodInferenceArray()
    : Method(),
      level_(0),
      array_size_(0),
      array_ghosts_(0),
      field_group_()
  { }

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodInferenceArray);

  /// Charm++ PUP::able migration constructor
  EnzoMethodInferenceArray (CkMigrateMessage *m)
    : Method (m),
      level_(0),
      array_size_(0),
      array_ghosts_(0),
      field_group_()
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) ;

  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "inference_array"; }

protected: // methods

  void compute_ (Block * block, enzo_float * Unew ) throw();

protected: // attributes

  /// Refinement level corresponding to the resolution of the inference array
  int level_;

  /// Size of the inference arrays excluding ghosts
  int array_size_;

  /// Depth of inference array ghost cells
  int array_ghosts_;

  /// Field group defining which fields to include in the inference
  /// array
  std::string field_group_;
};

#endif /* ENZO_ENZO_METHOD_INFERENCE_ARRAY_HPP */
