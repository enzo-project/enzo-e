// $Id: method_EnzoItBlocks.hpp 1942 2011-01-20 00:53:45Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef ENZO_ENZO_IT_BLOCKS_HPP
#define ENZO_ENZO_IT_BLOCKS_HPP

/// @file     method_EnzoItBlocks.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Implement operator++() to iterate over all Patches not just root
/// @todo     Implement iterators using templates to avoid void *
/// @todo     Move ItBlocks to Mesh component
/// @date     Tue Feb  1 16:46:01 PST 2011
/// @brief    [\ref Enzo] Declaration of the EnzoItBlocks iterator

class EnzoItBlocks : public Iterator {

  /// @class    EnzoItBlocks
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Iterator over Blocks in a Patch

public: // interface

  /// Create an EnzoItBlocks object
  EnzoItBlocks (Patch * patch, EnzoDescr * enzo) throw ();

  /// Delete the EnzoItBlocks object
  ~EnzoItBlocks () throw ();
  
  /// Iterate through all Patches in the patch (currently only root!)
  virtual void * operator++ ();

private: // attributes

  /// The Patch being iterated over
  Patch * patch_;

  /// Index to the current DataBlock
  int curr_;

  /// Enzo descriptor
  EnzoDescr * enzo_;
};

#endif /* ENZO_ENZO_IT_BLOCKS_HPP */
