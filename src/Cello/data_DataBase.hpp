// See LICENSE_CELLO file for license and copyright information

/// @file     data_DataBase.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-12-17
/// @brief    [\ref Data] Declaration of the DataBase class

#ifndef DATA_DATABASE_HPP
#define DATA_DATABASE_HPP

class DataBase {

  /// @class    DataBase
  /// @ingroup  Data
  /// @brief    [\ref Data] 

public: // interface

  /// Constructor
  DataBase() throw() {};

  /// Copy constructor
  DataBase(const DataBase & DataBase) throw() {};

  /// Assignment operator
  DataBase & operator= (const DataBase & DataBase) throw() {return *this;};

  /// Destructor
  virtual ~DataBase() throw() {};

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) {};
  
  //--------------------------------------------------

public: // virtual functions

  /// Return the number of bytes required to serialize the data object
  virtual int data_size () const = 0;

  /// Serialize the object into the provided empty memory buffer
  virtual char * save_data (char * buffer) const = 0;

  /// Restore the object from the provided initialized memory buffer data
  virtual char * load_data (char * buffer) = 0;

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

};

#endif /* DATA_DATABASE_HPP */

