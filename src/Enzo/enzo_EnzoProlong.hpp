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

  // /// Destructor
  // ~EnzoProlong() throw();

  // /// Copy constructor
  // EnzoProlong(const EnzoProlong & prolong) throw();

  // /// Assignment operator
  // EnzoProlong & operator= (const EnzoProlong & prolong) throw();

#ifdef CONFIG_USE_CHARM

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoProlong);

  /// CHARM++ migration constructor
  EnzoProlong(CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

#endif

  /// Prolong comm_block_Ht values to the child block given by (icx,icy,icz)
  virtual void apply (CommBlock        * comm_block_h, 
		      const CommBlock  * comm_block_Ht, 
		      const FieldDescr * field_descr,
		      int icx, int icy, int icz);

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Interpolation Method: see Enzo documenation
  int method_;

};

#endif /* FIELD_ENZO_PROLONG_HPP */

