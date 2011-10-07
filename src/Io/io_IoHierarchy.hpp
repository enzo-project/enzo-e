// See LICENSE_CELLO file for license and copyright information

#ifndef IO_IO_HIERARCHY_HPP
#define IO_IO_HIERARCHY_HPP

/// @file     io_IoHierarchy.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    [\ref Io] Declaration of the IoHierarchy class
///

class Hierarchy;

class IoHierarchy : public Io {

  /// @class    IoHierarchy
  /// @ingroup  Io
  /// @brief    [\ref Io] 

public: // interface

  /// Constructor
  IoHierarchy(const Hierarchy * hierarchy) throw();

  /// Return the ith metadata item associated with the Hierarchy object
  void meta_value 
  (int index, 
   void ** buffer, const char ** name, enum scalar_type * type,
   int * n0=0, int * n1=0, int * n2=0, int * n3=0, int * n4=0) throw();

  /// Return the ith data item associated with the Hierarchy object
  void data_value 
  (int index, 
   void ** buffer, const char ** name, enum scalar_type * type,
   int * n0=0, int * n1=0, int * n2=0, int * n3=0, int * n4=0) throw();

  
private: // functions

  const Hierarchy * hierarchy_;

private: // attributes


};

#endif /* IO_IO_HIERARCHY_HPP */

