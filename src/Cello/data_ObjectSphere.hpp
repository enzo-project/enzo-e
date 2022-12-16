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
  ObjectSphere(double center[3], double radius) throw()
    : Object(),
      radius_(radius)
  {
    center_[0] = center[0];
    center_[1] = center[1];
    center_[2] = center[2];
  };

  ObjectSphere() throw()
    : Object()
  { };

  /// Charm++ PUP::able declarations
  PUPable_decl(ObjectSphere);

  ObjectSphere (CkMigrateMessage *m)
    : Object(m),
      center_(),
      radius_(0)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;

    Object::pup(p);
    PUParray(p,center_,3);
    p | radius_;
  }

  ///--------------------
  /// PACKING / UNPACKING
  ///--------------------
  
  /// Return the number of bytes required to serialize the data object
  int data_size () const
  {
    //--------------------------------------------------
    //  1. determine buffer size (must be consistent with #3)
    //--------------------------------------------------

    int size = 0;

    SIZE_ARRAY_TYPE(size,double,center_,3);
    SIZE_SCALAR_TYPE(size,double,radius_);

    return size;
  }

  //----------------------------------------------------------------------

  /// Serialize the object into the provided empty memory buffer.
  /// Returns the next open position in the buffer to simplify
  /// serializing multiple objects in one buffer.
  char * save_data (char * buffer) const
  {
    char * pc = buffer;

    SAVE_ARRAY_TYPE(pc,double,center_,3);
    SAVE_SCALAR_TYPE(pc,double,radius_);

    ASSERT2 ("ObjectSphere::save_data()",
             "Expecting buffer size %d actual size %d",
             data_size(),(pc-buffer),
             (data_size() == (pc-buffer)));

    return pc;
  }

  //----------------------------------------------------------------------

  /// Restore the object from the provided initialized memory buffer data.
  /// Returns the next open position in the buffer to simplify
  /// serializing multiple objects in one buffer.
  char * load_data (char * buffer)
  {
    char * pc = buffer;

    LOAD_ARRAY_TYPE(pc,double,center_,3);
    LOAD_SCALAR_TYPE(pc,double,radius_);

    return pc;
  }

public: // virtual methods

  virtual void draw() { CkPrintf ("ObjectSphere::draw()\n"); }
  virtual void print(std::string msg) {
    CkPrintf ("ObjectSphere::print() %s center %g %g %g radius %g\n",
              msg.c_str(), center_[0],center_[1],center_[2],radius_); }

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

  double center_[3];
  double radius_;

};

#endif /* DATA_OBJECT_SPHERE_HPP */

