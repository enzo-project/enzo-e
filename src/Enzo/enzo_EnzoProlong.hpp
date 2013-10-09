// See LICENSE_CELLO file for license and copyright information

/// @file     field_EnzoProlong.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-09
/// @brief    [\ref Field] Declaration of the EnzoProlong class
///
/// This class serves to encapsulate Enzo's interpolate() function

#ifndef FIELD_ENZO_PROLONG_HPP
#define FIELD_ENZO_PROLONG_HPP

class EnzoProlong : public Prolong {

  /// @class    EnzoProlong
  /// @ingroup  Field
  /// @brief    [\ref ] 

public: // interface

  /// Constructor
  EnzoProlong(std::string method) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoProlong);

  /// CHARM++ migration constructor
  EnzoProlong(CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Prolong comm_block_Ht values to the child block given by (icx,icy,icz)

  virtual int apply 
  ( precision_type precision,
    void *       values_f, int nd3_f[3], int im3_f[3], int n3_f[3],
    const void * values_c, int nd3_c[3], int im3_c[3], int n3_c[3]);

private: // functions

  template <class T>
  int apply_
  ( T *       values_f, int nd3_f[3], int im3_f[3], int n3_f[3],
    const T * values_c, int nd3_c[3], int im3_c[3], int n3_c[3]);
  
private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Interpolation Method: see Enzo documenation
  int method_;

};

#endif /* FIELD_ENZO_PROLONG_HPP */

