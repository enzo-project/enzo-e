// See LICENSE_CELLO file for license and copyright information

/// @file     data_ObjectSphere.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-11-12
/// @brief    [\ref Data] Declaration of the ObjectSphere class

#ifndef DATA_OBJECT_SPHERE_HPP
#define DATA_OBJECT_SPHERE_HPP

class ObjectSphere : public Object {

  /// @class    ObjectSphere
  /// @ingroup  Data
  /// @brief    [\ref Data] 

public: // interface

  /// Constructor
  ObjectSphere() throw()
  { };

  /// Charm++ PUP::able declarations
  PUPable_decl(ObjectSphere);

  ObjectSphere (CkMigrateMessage *m)
    : Object(m),
      center_(),
      radius_(0)
  { }

  /// Copy constructor
  ObjectSphere(const ObjectSphere & ObjectSphere) throw()
  { };

  /// Assignment operator
  ObjectSphere & operator= (const ObjectSphere & ObjectSphere) throw()
  { };

  /// Destructor
  virtual ~ObjectSphere() throw()
  { };

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;

    Object::pup(p);
    PUParray(p,center_,3);
    p | radius_;
  }



public: // virtual methods

  virtual void draw() { CkPrintf ("ObjectSphere::draw()\n"); }
  virtual void print() { CkPrintf ("ObjectSphere::print()\n"); }

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

  double center_[3];
  double radius_;

};

#endif /* DATA_OBJECT_SPHERE_HPP */

