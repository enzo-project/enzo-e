// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_DataDescr.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 14:09:42 PDT 2010
/// @brief    [\ref Mesh] Declaration of the DataDescr class

#ifndef MESH_DATA_DESCR_HPP
#define MESH_DATA_DESCR_HPP

class FieldDescr;

class DataDescr {

  /// @class    DataDescr
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] Container class for all data descriptors (currently just fields)

public: // interface

  /// Initialize the DataDescr object
  DataDescr() throw();

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  ~DataDescr() throw();

  /// Copy constructor
  DataDescr(const DataDescr & data_descr) throw();

  /// Assignment operator
  DataDescr & operator= (const DataDescr & data_descr) throw();

  //----------------------------------------------------------------------

  /// Return the Field descriptor
  FieldDescr * field_descr () throw()
  { return field_descr_; };

  /// Return the constant Field descriptor
  const FieldDescr * field_descr () const throw()
  { return field_descr_; };

private: // attributes

  FieldDescr    * field_descr_;

};

#endif /* MESH_DATA_DESCR_HPP */

