// See LICENSE_CELLO file for license and copyright information

/// @file     io_ItFilePatch.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-03-27
/// @brief    [\ref Io] Declaration of the ItFilePatch class
///
/// This class is used to iterate over all Patches defined in
/// a data file

#ifndef IO_IT_FILE_PATCH_HPP
#define IO_IT_FILE_PATCH_HPP

class ItFilePatch {

  /// @class    ItFilePatch
  /// @ingroup  Io
  /// @brief    [\ref Io] 

public: // interface

  /// Create an ItFilePatch object
  ItFilePatch (const Input * input) throw ();

  /// Delete the ItFilePatch object
  virtual ~ItFilePatch () throw ();

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change
    p | *input_;
    p | index1_;
  }
#endif
  
  /// Iterate through all Patches in the File
  Patch * operator++ () throw();

  /// Return whether the iteration is complete
  bool done() const throw();

private: // attributes

  /// The Input object being iterated over
  Input * input_;

  /// Index of the current local Patch plus 1, or 0 if between iterations
  /// Always in the range 0 <= index1_ <= number of patches
  size_t index1_;

};

#endif /* IO_IT_FILE_PATCH_HPP */

