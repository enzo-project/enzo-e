// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MaskValue.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-03-28
/// @brief    [\ref Problem] Declaration of the MaskValue class
///

#ifndef PROBLEM_MASK_VALUE_HPP
#define PROBLEM_MASK_VALUE_HPP

class MaskValue : public Mask {

  /// @class    MaskValue
  /// @ingroup  Problem
  /// @brief    [\ref Problem] 

public: // interface

  /// Constructor
  MaskValue() throw() { };

  MaskValue(Parameters * parameters) throw();

  /// Destructor
  ~MaskValue() throw();

  /// Copy constructor
  MaskValue(const MaskValue & mask) throw();

  /// Assignment operator
  MaskValue & operator= (const MaskValue & mask) throw();

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
  }

  void evaluate
  (const Hierarchy * hierarchy,
   const CommBlock * comm_block,
   FieldBlock * field_block, 
   int index_field,  int index_value,
   std::string field_name, 
   const FieldDescr * field_descr,
   int n, bool * mask, bool * deflt,
   double * x, double * y, double * z, double t) throw () = 0;

  
private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change
  Parameters * parameters_;

  /// Masks for values: mask_[index_mask][i]
  bool ** mask_;

  /// Size of the masks
  int *nx_;
  int *ny_;
  int *nz_;

  /// number of fields
  int num_fields_;

  /// number of masked values per field
  int * num_masks_;

};

#endif /* PROBLEM_MASK_VALUE_HPP */

