// See LICENSE_CELLO file for license and copyright information

/// @file     field_EnzoRestrict.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-10
/// @brief    [\ref Field] Declaration of the EnzoRestrict class
///
/// This class serves to encapsulate Enzo's restriction operations

#ifndef FIELD_ENZO_RESTRICT_HPP
#define FIELD_ENZO_RESTRICT_HPP

class EnzoRestrict : public Restrict {

  /// @class    EnzoRestrict
  /// @ingroup  Field
  /// @brief    [\ref ] 

public: // interface

  /// Constructor
  EnzoRestrict(std::string method) throw();

#ifdef CONFIG_USE_CHARM

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoRestrict);

  /// CHARM++ migration constructor
  EnzoRestrict(CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

#endif

  /// Restrict comm_block_Ht values to the child block given by (icx,icy,icz)
  virtual void apply (CommBlock        * comm_block_c, 
		      const CommBlock  * comm_block_f, 
		      const FieldDescr * field_descr,
		      int icx, int icy, int icz);

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

};

#endif /* FIELD_ENZO_RESTRICT_HPP */

