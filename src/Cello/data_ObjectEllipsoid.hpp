// See LICENSE_CELLO file for license and copyright information

/// @file     data_ObjectEllipsoid.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-11-12
/// @brief    [\ref Data] Declaration of the ObjectEllipsoid class

#ifndef DATA_OBJECT_ELLIPSOID_HPP
#define DATA_OBJECT_ELLIPSOID_HPP

class ObjectEllipsoid : public Object {

  /// @class    ObjectEllipsoid
  /// @ingroup  Data
  /// @brief    [\ref Data] 

public: // interface

  /// Constructor
  ObjectEllipsoid() throw()
  { };

  /// Charm++ PUP::able declarations
  PUPable_decl(ObjectEllipsoid);

  ObjectEllipsoid (CkMigrateMessage *m)
    : Object(m),
      centers_(),
      radius_(0)
  { }

  /// Copy constructor
  ObjectEllipsoid(const ObjectEllipsoid & ObjectEllipsoid) throw()
  { };

  /// Assignment operator
  ObjectEllipsoid & operator= (const ObjectEllipsoid & ObjectEllipsoid) throw()
  { };

  /// Destructor
  virtual ~ObjectEllipsoid() throw()
  { };

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;

    Object::pup(p);
    PUParray(p,centers_[0],3);
    PUParray(p,centers_[1],3);
    p | radius_;
  }



public: // virtual methods

  virtual void draw() { CkPrintf ("ObjectEllipsoid::draw()\n"); }
  virtual void print() { CkPrintf ("ObjectEllipsoid::print()\n"); }

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

  double centers_[2][3];
  double radius_;

};

#endif /* DATA_OBJECT_ELLIPSOID_HPP */

