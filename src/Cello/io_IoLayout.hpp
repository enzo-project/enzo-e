// See LICENSE_CELLO file for license and copyright information

/// @file     io_IoLayout.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-03-28
/// @brief    [\ref Io] Declaration of the IoLayout class

#ifndef IO_IO_LAYOUT_HPP
#define IO_IO_LAYOUT_HPP

class Layout;

class IoLayout : public Io {

  /// @class    IoLayout
  /// @ingroup  Io
  /// @brief    [\ref Io] Class for interfacing between Layout and Output classes

public: // interface

  /// Constructor
  IoLayout(const Layout * layout) throw();

  /// Destructor
  virtual ~IoLayout() throw()
  {}

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {

    TRACEPUP;

    // NOTE: change this function whenever attributes change
    Io::pup(p);

    WARNING ("IoLayout::pup","skipping layout_");
    //    if (p.isUnpacking()) layout_ = new Layout;
    //    p | *layout_;
  }

#endif

#include "_io_Io_common.hpp"

  
private: // functions

  Layout * layout_;

private: // attributes


};

#endif /* IO_IO_LAYOUT_HPP */

