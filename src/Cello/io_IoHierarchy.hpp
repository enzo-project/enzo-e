// See LICENSE_CELLO file for license and copyright information

/// @file     io_IoHierarchy.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    [\ref Io] Declaration of the IoHierarchy class

#ifndef IO_IO_HIERARCHY_HPP
#define IO_IO_HIERARCHY_HPP

class Hierarchy;

class IoHierarchy : public Io {

  /// @class    IoHierarchy
  /// @ingroup  Io
  /// @brief    [\ref Io] Class for linking between Hierarchy and Output classes

public: // interface

  /// Constructor
  IoHierarchy(const Hierarchy * hierarchy) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(IoHierarchy);

  /// CHARM++ migration constructor
  IoHierarchy(CkMigrateMessage *m) : Io(m) {}

  /// Destructor
  virtual ~IoHierarchy () throw()
  {}

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {

    TRACEPUP;

    // NOTE: change this function whenever attributes change

    Io::pup(p);

    p | max_level_;
    PUParray(p,lower_,3);
    PUParray(p,upper_,3);
    PUParray(p,root_size_,3);
    PUParray(p,blocking_,3);
  }


  /// Return the ith metadata item associated with the object
  virtual void meta_value 
  (int index, 
   void ** buffer, std::string * name, int * type,
   int * nxd=0, int * nyd=0, int * nzd=0) throw();

  /// Copy the values to the object
  virtual void save_to (void *); 

private: // attributes

  int max_level_;
  double lower_[3];
  double upper_[3];
  int blocking_[3];
  int root_size_[3];

};

#endif /* IO_IO_HIERARCHY_HPP */

