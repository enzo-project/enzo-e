// See LICENSE_CELLO file for license and copyright information

/// @file     data_Object.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-11-12
/// @brief    [\ref Data] Declaration of the Object class

#ifndef DATA_OBJECT_HPP
#define DATA_OBJECT_HPP

class Object : public PUP::able {

  /// @class    Object
  /// @ingroup  Data
  /// @brief    [\ref Data] 

public: // interface

  /// Constructor
  Object() throw()
  { }

  /// Charm++ PUP::able declarations
  PUPable_abstract(Object);

  Object (CkMigrateMessage *m)
    : PUP::able(m)
  { }

  /// Copy constructor
  Object(const Object & Object) throw()
  { }

  /// Assignment operator
  Object & operator= (const Object & Object) throw()
  { }

  /// Destructor
  virtual ~Object() throw()
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;
    PUP::able::pup(p);
  }

public: // virtual methods

  virtual void draw() { }
  virtual void print(std::string) = 0;

  virtual int data_size () const = 0;
  virtual char * save_data (char * buffer) const = 0;
  virtual char * load_data (char * buffer) = 0;

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

};

#endif /* DATA_OBJECT_HPP */

