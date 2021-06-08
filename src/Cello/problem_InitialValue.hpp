// See LICENSE_CELLO file for license and copyright information

/// @file     problem_InitialValue.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Jan  4 19:26:38 PST 2011
/// @brief    [\ref Problem] Default initialization method

#ifndef PROBLEM_INITIAL_DEFAULT_HPP
#define PROBLEM_INITIAL_DEFAULT_HPP

class FieldData;

class InitialValue : public Initial {

  /// @class    InitialValue
  /// @ingroup  Method
  /// @brief    [\ref Method] Default initialization method

public: // interface

  /// CHARM++ constructor
  InitialValue() throw() { }

  /// Constructor
  InitialValue(Parameters * parameters,
               int cycle, double time) throw();

  PUPable_decl(InitialValue);

  InitialValue(CkMigrateMessage *m)
    : Initial (m),
      ic_pairs_()
  {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Read initialization values from Initial group in parameter file
  virtual void enforce_block (Block * block,
                              const Hierarchy * hierarchy
                              ) throw();

private: // functions

  void copy_values_ (FieldData * field_data,
		     double * value, int index_field,
		     int nx, int ny, int nz) throw();

  template<class T>
  void copy_precision_
  (T * field, int offset, double * value, int nx, int ny, int nz);

private: // attributes

  /// pairs of field names and Value objects
  std::vector<std::pair<std::string,Value>> ic_pairs_;

};

#endif /* PROBLEM_INITIAL_DEFAULT_HPP */
