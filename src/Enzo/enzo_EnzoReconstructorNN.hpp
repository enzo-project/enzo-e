// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoReconstructorNN.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Wed May 1 2019
/// @brief    [\ref Enzo] Implementation of Enzo's Nearest Neighbor
///           Reconstruction

#ifndef ENZO_ENZO_RECONSTRUCTOR_NN_HPP
#define ENZO_ENZO_RECONSTRUCTOR_NN_HPP

class EnzoReconstructorNN : public EnzoReconstructor
{
  /// @class    EnzoReconstructorNN
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates nearest neighbor reconstruction of
  ///           primitives at cell interfaces

public: // interface

  /// Create a new EnzoReconstructorNN
  EnzoReconstructorNN(std::vector<std::string> group_names)
    : EnzoReconstructor(group_names)
  { }

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoReconstructorNN);

  /// CHARM++ migration constructor for PUP::able
  EnzoReconstructorNN (CkMigrateMessage *m)
    : EnzoReconstructor(m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    EnzoReconstructor::pup(p);
  };

  void reconstruct_interface (Block *block, Grouping &prim_group,
			      Grouping &priml_group, Grouping &primr_group,
			      int dim, EnzoEquationOfState *eos,
			      int stale_depth);

  int total_staling_rate()
  { return 1; }

  int immediate_staling_rate()
  { return 0; }

};

#endif /* ENZO_ENZO_RECONSTRUCTOR_NN_HPP */
