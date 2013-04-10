// See LICENSE_CELLO file for license and copyright information

/// @file     io_IoPatch.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    [\ref Io] Declaration of the IoPatch class

#ifdef REMOVE_PATCH
#else /* REMOVE_PATCH */

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

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {

    TRACEPUP;

    // NOTE: change this function whenever attributes change

    Io::pup(p);

    WARNING ("IoPatch::pup","skipping patch_");
    //    p | *patch_;
    
  }
#endif

#include "_io_Io_common.hpp"

  
private: // attributes

  Patch * patch_;

};

#endif /* IO_IO_PATCH_HPP */

#endif /* REMOVE_PATCH */
