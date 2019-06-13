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
  /// @param reconstructable_groups A vector of reconstructable quantities
  ///     (that are also listed in FIELD_TABLE). These are used as group names
  ///     in the Grouping objects that store field names. In effect this is
  ///     used to register the quantities operated on by the Reconstructor
  /// @param passive_groups A vector with the names of the groups of passively
  ///     advected scalars that may be included. (If a group is listed here but
  ///     the Grouping object doesn't actually provide any fields in the group,
  ///     no problems are caused)
  /// @param solver The name of the Riemann solver to use. Valid names include
  ///     "nn" and "plm"
  static EnzoReconstructor* construct_reconstructor
    (std::vector<std::string> reconstructable_groups,
     std::vector<std::string> passive_groups, std::string solver);

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
  ///
  /// @param block holds data to be processed
  /// @param prim_group holds field names of the cell-centered reconstructable
  ///     primitives. This object is expected to have Grouping matching the
  ///     names registerred with the factory method
  /// @param priml_group,primr_group holds field names where the reconstructed
  ///     left/right face-centered primitives will be stored. The relevant
  ///     fields should be formally defined as cell-centered (to allow for
  ///     reuse along multiple dimensions). During the calculation, they are
  ///     treated as face-centered (without having values on the exterior faces
  ///     of the block). Consequentially there will be some unused space at the
  ///     end of the arrays.  
  /// @param dim Dimension along which to reconstruct interface values. Values
  ///     of 0, 1, and 2 correspond to the x, y, and z directions, respectively.
  /// @param eos Pointer to an instance of EnzoEquationOfState object
  /// @param stale_depth indicates the current stale_depth for the supplied
  ///     cell-centered quantities
  virtual void reconstruct_interface (Block *block, Grouping &prim_group,
				      Grouping &priml_group,
				      Grouping &primr_group, int dim,
				      EnzoEquationOfState *eos,
				      int stale_depth)=0;

  /// The rate amount by which the stale_depth increases after the current
  /// reconstructor is used to update the fluid over a (partial or full)
  /// time-step. If the fluid is update over a partial timestep before being
  /// updated, increases to the stale_depth are cummulative.
  ///
  /// stale_depth indicates the number of field entries from the outermost
  /// field value that the region including "stale" values (need to be
  /// refreshed) extends over. For cell-centered fields this is the number
  /// of cells away from the edge. Along the dimension(s) of face-centering,
  /// the region of valid (un-staled) entries for a face-centered field,
  /// with/without values on the exterior of the grid, will ALWAYS contain 1
  /// more/less valid field entry than the region of valid entries for a
  /// cell-centered field
  /// 
  /// the staling_rate can be decomposed into 2 parts: immediate and delayed
  /// - immediate staling rate is the amount by which stale_depth increases
  ///   immediately after performing reconstruction.
  /// - delayed staling rate is the amount by which stale_depth increases
  ///   after the fluid is updated over a (partial or full) time-step using the
  ///   fluxes computed from the reconstructed values.
  virtual int total_staling_rate()=0;

  /// immediate staling rate is the amount by which stale_depth increases
  /// immediately after performing reconstruction. This is equal to the number
  /// cells separating the first (last) cell interface that doesn't have both a
  /// left AND right reconstructed values from the edge of the grid. For
  /// interpolation like nearest neighbor, this is 0, while for interpolation
  /// like piecewise linear, this is 1.
  virtual int immediate_staling_rate()=0;

  /// delayed staling rate is the amount by which stale_depth increases after
  /// the fluid is updated over a (partial or full) time-step using the fluxes
  /// computed from the reconstructed values.
  int delayed_staling_rate()
  { return total_staling_rate() - immediate_staling_rate(); }

protected:
  /// list of the group_names which include all field names that will be
  /// reconstructed
  std::vector<std::string> group_names_;
};

#endif /* ENZO_ENZO_RECONSTRUCTOR_HPP */
