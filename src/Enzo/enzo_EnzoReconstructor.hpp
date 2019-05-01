// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoReconstructor.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Wed May 1 2019
/// @brief    [\ref Enzo] Implementation of Enzo's Reconstructor interface

#ifndef ENZO_ENZO_RECONSTRUCTOR_HPP
#define ENZO_ENZO_RECONSTRUCTOR_HPP

#include <pup_stl.h>

class EnzoReconstructor : public PUP::able
{
  /// @class    EnzoReconstructor
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates reconstruction of primitives at
  ///           cell interfaces

public: // interface

  /// Factory method for constructing EnzoReconstructor
  /// (The signature of this method may need to be modified)
  static EnzoReconstructor* construct_reconstructor
  (std::string name, const EnzoFieldConditions cond);

  /// Create a new EnzoReconstructor
  EnzoReconstructor(std::vector<std::string> group_names) throw()
    : group_names_(group_names)
  { }

  /// Virtual destructor
  virtual ~EnzoReconstructor()
  {  }

  /// CHARM++ PUP::able declaration
  PUPable_abstract(EnzoReconstructor);

  /// CHARM++ migration constructor for PUP::able
  EnzoReconstructor (CkMigrateMessage *m)
    : PUP::able(m)  
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    PUP::able::pup(p);
    p|group_names_;
  }

  /// Reconstructs the interface values
  /// @param block holds data to be processed
  /// @param prim_group holds field names of cell-centered primitives
  /// @param priml_group,primr_group holds field names where reconstructed
  ///  left/right face-centered primitives will be stored. The relevant fields
  ///  should be formally defined as cell-centered (to allow for reuse) and so
  ///  there will be some unused space at the end of the arrays.
  /// @param dim Dimension along which to reconstruct interface values. Values
  ///  of 0, 1, and 2 correspond to the x, y, and z directions, respectively.
  /// @param eos Instance of the fluid's EnzoEquationOfState object
  virtual void reconstruct_interface (Block *block, Grouping &prim_group,
				      Grouping &priml_group,
				      Grouping &primr_group, int dim,
				      EnzoEquationOfState *eos)=0;

  /// Returns the number of edge cells that become "stale" (need to be
  /// refreshed) each time the current reconstructor is used to update the fluid
  /// (if a fluid is first evolved over a half time-step and then a full
  /// time-step, cells become stale both times. The ghost depth must be at least
  /// as large as the total number of staled cells after the full time-step)
  virtual int get_staling_rate()=0;

protected:
  /// list of the group_names which include all field names that will be
  /// reconstructed
  std::vector<std::string> group_names_;
};

#endif /* ENZO_ENZO_RECONSTRUCTOR_HPP */
