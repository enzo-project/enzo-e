// See LICENSE_CELLO file for license and copyright information

/// @file     field_EnzoProlongMC1.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-09
/// @brief    [\ref Field] Declaration of the EnzoProlongMC1 class
///
/// This class serves to encapsulate Enzo's interpolate() function

#ifndef FIELD_ENZO_PROLONG_MC1_HPP
#define FIELD_ENZO_PROLONG_MC1_HPP

class EnzoProlongMC1 : public Prolong {

  /// @class    EnzoProlongMC1
  /// @ingroup  Field
  /// @brief    [\ref ] 

public: // interface

  /// Constructor
  EnzoProlongMC1(std::string method) throw();

#ifdef CONFIG_USE_CHARM

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoProlongMC1);

  /// CHARM++ migration constructor
  EnzoProlongMC1(CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

#endif

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

#endif /* FIELD_ENZO_PROLONG_MC1_HPP */

