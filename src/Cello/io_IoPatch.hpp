// See LICENSE_CELLO file for license and copyright information

/// @file     io_IoPatch.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    [\ref Io] Declaration of the IoPatch class

#ifndef IO_IO_PATCH_HPP
#define IO_IO_PATCH_HPP

class Patch;

class IoPatch : public Io {

  /// @class    IoPatch
  /// @ingroup  Io
  /// @brief    [\ref Io] Class for interfacing between Patch and Output classes

public: // interface

  /// Constructor
  IoPatch(const Patch * patch) throw();

  /// Destructor
  virtual ~IoPatch() throw()
  {}

  /// Return the ith metadata item associated with the Patch object
  void meta_value 
  (int index, 
   void ** buffer, std::string * name, enum scalar_type * type,
   int * n0=0, int * n1=0, int * n2=0, int * n3=0, int * n4=0) throw();

  /// Return the ith data item associated with the Patch object
  void data_value 
  (int index, 
   void ** buffer, std::string * name, enum scalar_type * type,
   int * n0=0, int * n1=0, int * n2=0, int * n3=0, int * n4=0) throw();

  
private: // functions

  const Patch * patch_;

private: // attributes


};

#endif /* IO_IO_PATCH_HPP */

